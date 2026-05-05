"""Module 3: Population junction quantification.

Quantifies known tandem duplication junctions in population sequencing data.
Uses exact substring matching (grep -F) to count reads containing the junction,
then calculates Junction Read Ratio (JRR) proportional to duplication frequency.

JRR interpretation:
  JRR ≈ 0.0   → junction not present (no duplicated cells)
  JRR ≈ 0.5   → ~50% of cells carry the duplication
  JRR ≈ 1.0   → all cells carry the duplication

The formula normalizes junction-spanning reads by the expected number of
reads at 100% duplication, given the junction length and WT coverage:
  expected_at_100pct = WT_coverage × (read_length - junction_length + 1) / read_length
  JRR = spanning_reads / expected_at_100pct

Usage: tandem -r reference.fna -i R1.fq -I R2.fq -pop -m metadata.tsv
"""

import logging
import os
import json
import subprocess
import time
from pathlib import Path

from . import utils
from . import alignment
from .junction import (
    KmerIndex,
    build_junction_reference,
    reverse_complement,
    extend_to_unique,
    classify_junction,
    classify_single_copy_junction,
)

logger = logging.getLogger("tandem")

# The output unit name
JRR_UNIT = "JRR"  # Junction Read Ratio


def _calculate_mean_coverage(bam_path, sequences):
    """Calculate mean coverage across the entire reference genome.

    Uses samtools depth on the sorted/indexed BAM file. Averages the
    depth across all positions in all sequences in the reference.

    Args:
        bam_path: path to sorted, indexed BAM file
        sequences: dict of seq_id -> sequence string (to get genome length)

    Returns:
        mean coverage (float)
    """
    result = utils.run_command(
        ["samtools", "depth", "-a", str(bam_path)],
        description="Reading whole-genome coverage from BAM",
        check=False,  # may produce many lines but not fail
    )

    total_depth = 0
    n_positions = 0
    if result.stdout:
        for line in result.stdout.strip().split('\n'):
            if not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            try:
                depth = int(parts[2])
                total_depth += depth
                n_positions += 1
            except (ValueError, IndexError):
                continue

    if n_positions == 0:
        logger.warning("  No coverage data found in BAM — returning 0")
        return 0.0

    return total_depth / n_positions


def _detect_read_length(fastq_path, n_reads=1000):
    """Detect read length by sampling first N reads from FASTQ."""
    lengths = []
    if str(fastq_path).endswith('.gz'):
        cmd = f"zcat {fastq_path} | head -n {n_reads * 4}"
    else:
        cmd = f"head -n {n_reads * 4} {fastq_path}"

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    for i, line in enumerate(result.stdout.split('\n')):
        if i % 4 == 1 and line.strip():
            lengths.append(len(line.strip()))
        if len(lengths) >= n_reads:
            break

    if not lengths:
        return 150  # default
    # Use median
    lengths.sort()
    return lengths[len(lengths) // 2]


def _count_exact_matches(junction_seqs, read1, read2, output_dir):
    """Count exact matches for each junction sequence in reads.

    Uses grep -F (Aho-Corasick) for fast multi-pattern exact matching.
    Searches both forward and reverse-complement of each junction.

    Returns:
        dict of {junction_sequence: count}
    """
    t0 = time.time()

    # Write junction sequences + reverse complements to patterns file
    patterns_file = os.path.join(output_dir, "junction_patterns.txt")
    seq_to_original = {}  # map RC back to original

    with open(patterns_file, 'w') as f:
        for seq in junction_seqs:
            f.write(seq + '\n')
            seq_to_original[seq] = seq
            rc = reverse_complement(seq)
            if rc != seq:
                f.write(rc + '\n')
                seq_to_original[rc] = seq

    # Extract read sequences from FASTQ to temp file
    reads_file = os.path.join(output_dir, "read_seqs.txt")
    for i, read_path in enumerate([read1, read2]):
        if read_path is None:
            continue
        mode = 'w' if i == 0 else 'a'
        if str(read_path).endswith('.gz'):
            cmd = f"zcat {read_path} | awk 'NR%4==2'"
        else:
            cmd = f"awk 'NR%4==2' {read_path}"

        logger.info(f"    Extracting read sequences ({i+1})")
        with open(reads_file, mode) as outf:
            subprocess.run(cmd, shell=True, stdout=outf, check=True)

    logger.info(f"    Prepared patterns and reads in {time.time()-t0:.1f}s")

    # Run grep -oFf: exact multi-pattern matching
    t0 = time.time()
    grep_result = utils.run_command(
        ['grep', '-oFf', patterns_file, reads_file],
        description="Exact match search (grep -F)",
        check=False,
    )

    # Count per-junction matches
    read_counts = {}
    if grep_result.stdout and grep_result.stdout.strip():
        for match in grep_result.stdout.strip().split('\n'):
            match = match.strip()
            if not match:
                continue
            original = seq_to_original.get(match, match)
            read_counts[original] = read_counts.get(original, 0) + 1

    elapsed = time.time() - t0
    total_hits = sum(read_counts.values())
    logger.info(
        f"    Exact match done: {len(read_counts)} junctions with hits, "
        f"{total_hits} total reads, {elapsed:.1f}s"
    )

    # Clean up temp files
    for f in [patterns_file, reads_file]:
        try:
            os.remove(f)
        except OSError:
            pass

    return read_counts


def _generate_candidates_worker(args):
    """Worker: test a chunk of (name, s, e) positions.

    Each worker builds its own KmerIndex to avoid pickling issues.
    Returns list of (name, s, e, seq, k) for valid unique junctions.
    """
    sequences, tasks, read_length = args

    # Build KmerIndex in this worker process
    worker_index = KmerIndex(sequences)
    # Use the largest sequence as main
    main_seq_id = max(sequences, key=lambda k: len(sequences[k]))
    main_seq = sequences[main_seq_id]

    results = []
    for name, s, e in tasks:
        result = extend_to_unique(main_seq, s, e, worker_index)
        if result is None:
            continue
        if len(result["sequence"]) > read_length:
            continue
        results.append((name, s, e, result["sequence"], result["k"]))
    return results


def run_module3(ref_fasta, read1, read2=None, metadata_path=None,
                output_dir=".", threads=1, flank=25,
                hr_complexity_filter=True):
    """Run the complete module 3 pipeline.

    For each metadata entry, searches a small neighborhood (±flank bp)
    around the provided coordinates to find the actual junction. This
    handles the typical 1-20bp uncertainty in reported coordinates from
    microhomology, different tools, or coordinate conventions.

    Step 1: Parse metadata to get approximate junction coordinates.
    Step 2: For each entry, generate candidate junctions in ±flank window.
    Step 3: Count exact matches in reads (grep -F).
    Step 4: Pick top-hit junction per metadata entry.
    Step 5: Calculate WT coverage and JRR; classify HR vs non-HR.

    Args:
        ref_fasta: path to reference FASTA
        read1: path to forward reads
        read2: path to reverse reads (optional)
        metadata_path: path to metadata TSV file (required)
        output_dir: output directory
        threads: number of threads (for WT mapping)
        flank: ±flank bp search window around provided coordinates (default 25)

    Returns:
        list of junction quantification result dicts
    """
    logger.info("=" * 60)
    logger.info("Module 3: Population junction quantification")
    logger.info("=" * 60)

    output_dir = utils.ensure_dir(output_dir)

    if metadata_path is None:
        logger.error(
            "Module 3 requires a metadata file (-m). "
            "Run module 2 (-iso) first to discover junctions, "
            "or provide metadata manually."
        )
        return []

    # Load reference
    logger.info(f"Loading reference: {ref_fasta}")
    sequences, headers = utils.load_fasta(ref_fasta)
    main_seq_id = max(sequences, key=lambda k: len(sequences[k]))
    main_seq = sequences[main_seq_id]
    logger.info(f"  Primary sequence: {main_seq_id} ({len(main_seq):,} bp)")

    # Build k-mer index
    logger.info("Building k-mer index...")
    kmer_index = KmerIndex(sequences)

    # Detect read length
    read_length = _detect_read_length(read1)
    logger.info(f"  Detected read length: {read_length}bp")

    # Step 1: Parse metadata
    logger.info(f"Step 1: Parsing metadata from {metadata_path}")
    entries = utils.parse_metadata(metadata_path)

    if not entries:
        logger.error("No valid entries in metadata file.")
        return []

    # Step 2: Build junction candidates within flank window (parallelized)
    logger.info(f"Step 2: Building junction candidates (flank=±{flank}bp)")

    # Collect all (entry_name, s, e) combinations to test in parallel
    from concurrent.futures import ProcessPoolExecutor, as_completed
    import time as _time

    tasks = []  # list of (entry_name, s, e)
    entry_info = {}  # name -> (start0, end0, provided_seq)
    entry_to_candidates = {}
    seq_to_entries = {}
    failed_junctions = []

    for entry in entries:
        name = entry["name"]
        start0 = entry["start"]
        end0 = entry["end"]
        provided_seq = entry.get("sequence")
        entry_info[name] = (start0, end0, provided_seq)
        entry_to_candidates[name] = []

        # Handle user-provided sequence directly
        if provided_seq and kmer_index.is_absent(provided_seq) and len(provided_seq) <= read_length:
            entry_to_candidates[name].append(
                (start0, end0, provided_seq, len(provided_seq) // 2)
            )
            continue

        # Collect tasks for this entry
        for ds in range(-flank, flank + 1):
            for de in range(-flank, flank + 1):
                s = start0 + ds
                e = end0 + de
                if e <= s or (e - s) < 100:
                    continue
                tasks.append((name, s, e))

    total_tasks = len(tasks)
    logger.info(
        f"  Testing {total_tasks:,} candidate positions across "
        f"{len(entries)} entries using {threads} processes"
    )

    if tasks:
        t0 = _time.time()
        # Split tasks into chunks (10× more chunks than threads for load balancing)
        n_chunks = min(threads * 10, len(tasks))
        chunks = []
        for i in range(n_chunks):
            chunk_start = i * len(tasks) // n_chunks
            chunk_end = (i + 1) * len(tasks) // n_chunks
            chunks.append((sequences, tasks[chunk_start:chunk_end], read_length))

        all_results = []  # list of (name, s, e, seq, k)
        completed = 0
        log_interval = max(1, len(chunks) // 10)

        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {
                executor.submit(_generate_candidates_worker, chunk): i
                for i, chunk in enumerate(chunks)
            }
            for future in as_completed(futures):
                all_results.extend(future.result())
                completed += 1
                if completed % log_interval == 0 or completed == len(chunks):
                    elapsed = _time.time() - t0
                    pct = 100 * completed / len(chunks)
                    logger.info(
                        f"    Progress: {completed}/{len(chunks)} chunks "
                        f"({pct:.0f}%), {len(all_results)} valid junctions, "
                        f"{elapsed:.1f}s"
                    )

        # Assign results to entries
        for name, s, e, seq, k in all_results:
            entry_to_candidates[name].append((s, e, seq, k))

    # Deduplicate per entry and build seq_to_entries
    for name, candidates in list(entry_to_candidates.items()):
        if not candidates:
            start0, end0, _ = entry_info[name]
            original = next(en for en in entries if en["name"] == name)
            logger.warning(
                f"  {name}: No unique junction found within ±{flank}bp. Skipping."
            )
            failed_junctions.append({
                "name": name,
                "start": original["user_start_1based"],
                "end": original["user_end_1based"],
                "reason": "no_unique_junction_in_flank",
            })
            del entry_to_candidates[name]
            continue

        seen = set()
        unique_candidates = []
        for c in candidates:
            if c[2] not in seen:
                seen.add(c[2])
                unique_candidates.append(c)
        entry_to_candidates[name] = unique_candidates

        for s, e, seq, k in unique_candidates:
            if seq not in seq_to_entries:
                seq_to_entries[seq] = []
            seq_to_entries[seq].append((name, s, e, k))

        logger.info(
            f"  {name}: {len(unique_candidates)} unique candidate junctions"
        )

    if not entry_to_candidates:
        logger.error("No junction candidates could be built for any entry.")
        _save_results([], failed_junctions, output_dir)
        return []

    # Step 3: Count exact matches (grep -F)
    logger.info(f"Step 3: Counting junction reads by exact match "
                f"({len(seq_to_entries)} unique sequences)")
    read_counts = _count_exact_matches(
        list(seq_to_entries.keys()), read1, read2, output_dir
    )

    # Step 4: Pick top-hit junction per entry
    logger.info("Step 4: Selecting top-hit junction per entry")
    entry_to_best = {}
    for name, candidates in entry_to_candidates.items():
        best = None
        best_count = -1
        for s, e, seq, k in candidates:
            count = read_counts.get(seq, 0)
            if count > best_count:
                best_count = count
                best = (s, e, seq, k, count)
        entry_to_best[name] = best

    # Step 5: Calculate WT mean coverage by mapping reads to reference
    # We use the whole-genome average coverage (from mapped reads) as
    # the baseline. This represents the "effective sequencing depth" —
    # only reads that belong to this genome are counted.
    logger.info("Step 5: Calculating WT mean coverage (mapping reads to reference)")
    wt_bam = os.path.join(output_dir, "reads_to_wt.sorted.bam")
    # Use conservative thread count for mapping to avoid memory/pipe issues
    map_threads = min(threads, 16)
    alignment.map_reads_to_reference(
        ref_fasta, read1, read2, output_bam=wt_bam, threads=map_threads
    )

    # Calculate mean coverage across whole genome using samtools depth
    logger.info("  Calculating genome-wide mean coverage")
    wt_cov_genome = _calculate_mean_coverage(wt_bam, sequences)
    logger.info(f"  Mean WT coverage: {wt_cov_genome:.2f}x")

    # Step 6: Calculate JRR and classify mechanism
    logger.info("Step 6: Calculating JRR and classifying mechanism")
    results = []
    for name, best in entry_to_best.items():
        if best is None:
            continue
        s, e, seq, k, spanning = best
        junc_len = len(seq)

        # Retrieve original user-provided coordinates
        original = next(en for en in entries if en["name"] == name)
        orig_start = original["start"]            # 0-based internal
        orig_end = original["end"]                # 1-based inclusive internal
        user_start_1based = original["user_start_1based"]  # for display
        user_end_1based = original["user_end_1based"]

        # Use whole-genome mean coverage as baseline
        wt_cov = wt_cov_genome

        # JRR calculation
        if wt_cov > 0 and read_length > junc_len:
            expected_at_100pct = wt_cov * (read_length - junc_len + 1) / read_length
            jrr = spanning / expected_at_100pct if expected_at_100pct > 0 else 0.0
        else:
            expected_at_100pct = 0.0
            jrr = 0.0

        # Classify junction using single-copy reference analysis.
        # The reference has only one copy of the duplicated region;
        # classify_single_copy_junction correctly measures MH (end-of-Y
        # vs start-of-Y) and HR (R-Y-R flanking structure).
        classification = classify_single_copy_junction(
            main_seq, s, e, hr_complexity_filter=hr_complexity_filter,
        )

        # Convert internal 0-based start / Python exclusive end to 1-based inclusive for output
        display_start = s + 1
        display_end = e  # 1-based inclusive == Python exclusive

        result = {
            "name": name,
            "seq_id": main_seq_id,
            "provided_start": user_start_1based,
            "provided_end": user_end_1based,
            "dup_start": display_start,
            "dup_end": display_end,
            "dup_size": e - s,
            "position_shift_start": display_start - user_start_1based,
            "position_shift_end": display_end - user_end_1based,
            "junction_id": f"{name}_s{display_start}_e{display_end}",
            "junction_k": k,
            "junction_length": junc_len,
            "spanning_reads": spanning,
            "wt_mean_coverage": round(wt_cov, 2),
            "expected_reads_100pct": round(expected_at_100pct, 2),
            "read_length": read_length,
            JRR_UNIT: round(jrr, 4),
            "is_hr_signature": classification["is_hr_signature"],
            "hr_match_len": classification["hr_match_len"],
            "hr_identity": classification["hr_identity"],
            "microhomology_bp": classification["microhomology_bp"],
            "microhomology_seq": classification["microhomology_seq"],
        }

        results.append(result)

        mech = "HR" if classification["is_hr_signature"] else f"non-HR (MH={classification['microhomology_bp']}bp)"
        shift = ""
        if display_start != user_start_1based or display_end != user_end_1based:
            shift = (f", shifted ({display_start - user_start_1based:+d}, "
                     f"{display_end - user_end_1based:+d})")
        logger.info(
            f"  {name}: {display_start:,}-{display_end:,} {mech}{shift}, "
            f"{spanning} reads, WT cov={wt_cov:.1f}x, "
            f"{JRR_UNIT}={jrr:.4f}"
        )

    # Save results
    _save_results(results, failed_junctions, output_dir)

    return results


def _save_results(results, failed_junctions, output_dir):
    """Save module 3 results."""
    # TSV
    tsv_path = os.path.join(output_dir, "junction_quantification.tsv")
    with open(tsv_path, "w") as f:
        f.write(f"# Tandem Module 3: Population junction quantification\n")
        f.write(f"# {JRR_UNIT} = Junction Read Ratio (~0.0 = absent, "
                f"~1.0 = all cells duplicated)\n")
        f.write(f"# Formula: JRR = spanning_reads / (wt_coverage × "
                f"(read_length - junction_length + 1) / read_length)\n")

        if results:
            header = list(results[0].keys())
            f.write("\t".join(header) + "\n")
            for r in results:
                row = [str(r.get(h, "")) for h in header]
                f.write("\t".join(row) + "\n")
        else:
            f.write("# No junctions quantified\n")

    logger.info(f"  Results saved to {tsv_path}")

    # JSON
    json_path = os.path.join(output_dir, "junction_quantification.json")
    output = {
        "unit": JRR_UNIT,
        "unit_description": (
            "Junction Read Ratio - proportional to duplication frequency in "
            "the population. ~0.0 means junction absent, ~1.0 means all cells "
            "carry the duplication. Calculated as observed spanning reads "
            "divided by expected spanning reads at 100% duplication."
        ),
        "results": results,
        "failed_junctions": failed_junctions,
    }
    with open(json_path, "w") as f:
        json.dump(output, f, indent=2)

    logger.info(f"  JSON results saved to {json_path}")
