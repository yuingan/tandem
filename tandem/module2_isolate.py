"""Module 2: Isolate resequencing junction discovery.

Discovers tandem duplication junctions from whole-genome resequencing data
of an isolate against a known reference genome.

Usage: tandem -r reference.fna -i R1.fq -I R2.fq -iso
"""

import logging
import os
import json
from pathlib import Path

import numpy as np

from . import utils
from . import alignment
from . import coverage as cov_module
from .junction import (
    KmerIndex,
    generate_all_junction_candidates,
    classify_junction,
    classify_single_copy_junction,
    extend_to_unique,
)

logger = logging.getLogger("tandem")


def run_module2(ref_fasta, read1, read2=None, output_dir=".",
                threads=1, window=200, flank=200,
                junction_window=None,
                manual_start=None, manual_end=None,
                max_k=75, hr_complexity_filter=True):
    """Run Module 2: Isolate junction confirmation.

    Workflow:
      1. Map reads to reference, calculate coverage.
      2. If -s/-e coordinates provided: generate candidate junctions
         from all start × end combinations, confirm by exact match
         in reads, classify mechanism.
      3. If no -s/-e: generate interactive coverage plot and exit.
         User examines the plot to identify candidate duplication
         boundaries, then re-runs with -s/-e.

    Args:
        ref_fasta: path to reference FASTA
        read1: path to forward reads
        read2: path to reverse reads (optional)
        output_dir: output directory
        threads: number of threads
        window: sliding window size for coverage calculation
        flank: flank size for junction candidate generation
        junction_window: search window around each coordinate.
            If 0 (--precise), search exact position only.
            If None, defaults to `window`.
        manual_start: list of candidate start positions (0-based internal)
        manual_end: list of candidate end positions (0-based internal)
        max_k: maximum junction extension length
        hr_complexity_filter: reject HR calls from low-complexity windows

    Returns:
        list of confirmed junction dicts
    """
    logger.info("=" * 60)
    logger.info("Module 2: Isolate junction confirmation")
    logger.info("=" * 60)

    output_dir = utils.ensure_dir(output_dir)

    # Junction search window: 0 with --precise, otherwise same as coverage window
    jw = junction_window if junction_window is not None else window

    # Load reference
    logger.info(f"Loading reference: {ref_fasta}")
    sequences, headers = utils.load_fasta(ref_fasta)
    total_bp = sum(len(s) for s in sequences.values())
    logger.info(f"  {len(sequences)} sequence(s), {total_bp:,} bp total")

    # Build k-mer index
    logger.info("Building k-mer index...")
    kmer_index = KmerIndex(sequences)

    # Step 1: Map reads to reference and calculate coverage
    logger.info("Step 1: Mapping reads to reference")
    bam_path = os.path.join(output_dir, "reads_to_ref.sorted.bam")
    alignment.map_reads_to_reference(
        ref_fasta, read1, read2, output_bam=bam_path, threads=threads
    )

    logger.info("Step 2: Calculating coverage")
    all_coverage = {}
    all_medians = {}
    for seq_id in sequences:
        logger.info(f"  Calculating coverage for {seq_id}...")
        cov_data = cov_module.calculate_coverage_from_bam(bam_path, seq_id)
        if seq_id in cov_data and len(cov_data[seq_id]) > 0:
            all_coverage[seq_id] = cov_data[seq_id]
            all_medians[seq_id] = float(np.median(cov_data[seq_id]))

    # --- No coordinates provided: coverage plot only ---
    if manual_start is None or manual_end is None:
        logger.info("")
        logger.info("No candidate coordinates provided (-s/-e).")
        logger.info("Generating coverage plot for manual inspection.")
        logger.info("")
        logger.info("Workflow:")
        logger.info("  1. Examine the HTML coverage plot below")
        logger.info("  2. Identify candidate duplication start/end boundaries")
        logger.info("  3. Re-run with: tandem -r REF -i R1 -I R2 -iso \\")
        logger.info("       -s <start1> <start2> ... -e <end1> <end2> ... \\")
        logger.info("       --precise -flank 400 -t THREADS -o OUTPUT")

        for seq_id in sequences:
            if seq_id not in all_coverage:
                continue
            genome_median = all_medians.get(seq_id, 0)
            cov_module.generate_coverage_plot(
                all_coverage[seq_id],
                seq_id=seq_id,
                output_dir=output_dir,
                genome_median=genome_median,
                elevated_regions=[],
            )
            logger.info(f"  Coverage plot saved to {output_dir}/")

        return []

    # --- Coordinates provided: junction confirmation ---
    starts = manual_start if isinstance(manual_start, list) else [manual_start]
    ends = manual_end if isinstance(manual_end, list) else [manual_end]

    target_seq_id = None
    for seq_id, seq in sequences.items():
        if max(ends) <= len(seq):
            target_seq_id = seq_id
            break
    if target_seq_id is None:
        target_seq_id = list(sequences.keys())[0]

    # Generate all start × end combinations
    regions = []
    for s in starts:
        for e in ends:
            if e > s:
                regions.append({
                    "start": s,
                    "end": e,
                    "mean_coverage": 0,
                    "fold_change": 0,
                    "size": e - s,
                    "source": "manual",
                })

    elevated_regions = {target_seq_id: regions}

    logger.info(
        f"  Candidate coordinates: {len(starts)} start(s) × {len(ends)} end(s) "
        f"= {len(regions)} combinations"
    )
    for r in regions:
        logger.info(
            f"    {target_seq_id}:{r['start']+1:,}-{r['end']:,} "
            f"({r['size']:,} bp)"
        )

    # Generate coverage plot with candidate regions marked
    if target_seq_id in all_coverage:
        genome_median = all_medians.get(target_seq_id, 0)
        cov_module.generate_coverage_plot(
            all_coverage[target_seq_id],
            seq_id=target_seq_id,
            output_dir=output_dir,
            genome_median=genome_median,
            elevated_regions=regions,
        )
    # Step 3-5: Generate and confirm junction candidates
    logger.info("Step 3: Generating junction candidates from boundary pairs")
    all_confirmed = []

    ref_seq = sequences[target_seq_id]

    # Generate boundary pairs (all start × end combinations, end > start)
    pairs = []
    for s in starts:
        for e in ends:
            if e > s and (e - s) >= 100:
                pairs.append((s, e))

    logger.info(f"  Testing {len(pairs)} boundary combinations")

    # Generate candidates for each pair
    all_pair_candidates = []
    for pair_idx, (s, e) in enumerate(pairs):
        logger.info(
            f"\n  Pair {pair_idx + 1}/{len(pairs)}: "
            f"start={s+1:,}, end={e:,} ({e - s:,} bp)"
        )

        candidates = generate_all_junction_candidates(
            ref_seq,
            s_approx=s,
            e_approx=e,
            window=jw,
            flank=flank,
            kmer_index=kmer_index,
            threads=threads,
            sequences=sequences,
            max_k=max_k,
        )
        all_pair_candidates.extend(candidates)

    # Deduplicate across all pairs
    seen = set()
    deduped = []
    for cand in all_pair_candidates:
        if cand["sequence"] not in seen:
            seen.add(cand["sequence"])
            deduped.append(cand)

    logger.info(f"\n  Total unique candidates: {len(deduped)}")

    if deduped:
        confirmed = _process_candidates(
            deduped, ref_seq, target_seq_id, kmer_index,
            read1, read2, output_dir, 0, threads,
            hr_complexity_filter=hr_complexity_filter,
        )
        all_confirmed.extend(confirmed)

    # Save results
    logger.info(f"\n  Total confirmed junctions: {len(all_confirmed)}")
    _save_results(all_confirmed, output_dir)

    return all_confirmed


def _process_candidates(candidates, ref_seq, seq_id, kmer_index,
                        read1, read2, output_dir, group_idx, threads,
                        hr_complexity_filter=True):
    """Confirm junction candidates by exact substring match in raw reads.

    Uses grep -F (Aho-Corasick) for fast multi-pattern exact matching.
    Each junction sequence must appear as an exact continuous substring
    within a read. No approximate matching, no soft-clipping.

    grep gives accurate counts with no multi-mapping confusion because
    each junction is a unique sequence — a read either contains it or not.

    Returns:
        list of confirmed junction dicts
    """
    if not candidates:
        logger.info("    No valid junction candidates found")
        return []

    import time
    import subprocess
    from tandem.junction import reverse_complement

    # Collect unique junction sequences
    junction_seqs = {}  # seq -> list of candidates
    for cand in candidates:
        # junction_id uses 1-based coordinates for consistency with output
        junc_id = (
            f"junc_{len(junction_seqs):04d}"
            f"_s{cand['start']+1}_e{cand['end']}_k{cand['k']}"
        )
        cand["junction_id"] = junc_id
        seq = cand["sequence"]
        if seq not in junction_seqs:
            junction_seqs[seq] = []
        junction_seqs[seq].append(cand)

    logger.info(
        f"    Confirming {len(junction_seqs)} unique junction sequences "
        f"by exact match in reads (grep -F)"
    )

    t0 = time.time()

    # Write junction sequences + reverse complements to patterns file
    patterns_file = os.path.join(output_dir, f"junction_patterns_{group_idx}.txt")
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
    reads_file = os.path.join(output_dir, f"read_seqs_{group_idx}.txt")
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
    # -o: output only matching part (so we know which pattern matched)
    # -F: fixed strings (exact match, no regex)
    # -f: patterns from file
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

    if not read_counts:
        logger.info(
            f"    No exact matches found for any junction candidate. "
            f"Region is likely noise — skipping."
        )
        return []

    # Build confirmed junctions
    confirmed = []
    for seq, count in sorted(read_counts.items(), key=lambda x: -x[1]):
        cand_list = junction_seqs.get(seq, [])
        if not cand_list:
            continue
        cand = cand_list[0]

        # Classify the junction using single-copy reference analysis.
        # The reference has only one copy; classify_single_copy_junction
        # correctly measures MH (end-of-Y vs start-of-Y) and HR (R-Y-R
        # structure flanking the single copy).
        classification = classify_single_copy_junction(
            ref_seq, cand["start"], cand["end"],
            hr_complexity_filter=hr_complexity_filter,
        )

        # Convert internal 0-based start / Python exclusive end to 1-based inclusive for output
        display_start = cand["start"] + 1
        display_end = cand["end"]

        result = {
            "seq_id": seq_id,
            "dup_start": display_start,
            "dup_end": display_end,
            "dup_size": cand.get("dup_size", cand["end"] - cand["start"]),
            "junction_id": cand["junction_id"],
            "junction_sequence": cand["sequence"],
            "junction_k": cand["k"],
            "total_reads": count,
            "hq_reads": count,
            "spanning_reads": count,
            "is_hr_signature": classification["is_hr_signature"],
            "hr_match_len": classification["hr_match_len"],
            "hr_identity": classification["hr_identity"],
            "hr_scenario": classification["hr_scenario"],
            "microhomology_bp": classification["microhomology_bp"],
            "microhomology_seq": classification["microhomology_seq"],
        }

        confirmed.append(result)

    # Sort by spanning reads
    confirmed.sort(key=lambda x: x["spanning_reads"], reverse=True)

    # Log confirmed junctions
    for hit in confirmed:
        hr = hit.get("is_hr_signature", False)
        mh = hit.get("microhomology_bp", 0)
        logger.info(
            f"    CONFIRMED: {hit['seq_id']}:{hit['dup_start']:,}-{hit['dup_end']:,} "
            f"({hit['dup_size']:,} bp) "
            f"{'HR' if hr else f'MH={mh}bp'} "
            f"[{hit['spanning_reads']} spanning reads, k={hit['junction_k']}]"
        )

    return confirmed


def _filter_confirmed_junctions(candidates, read_counts, ref_seq, seq_id,
                                 kmer_index, min_spanning=1,
                                 hr_complexity_filter=True):
    """Filter junction candidates by read support and classify.

    Args:
        candidates: list of junction candidate dicts
        read_counts: dict from count_junction_reads
        ref_seq: reference sequence
        seq_id: sequence ID
        kmer_index: KmerIndex object
        min_spanning: minimum spanning reads to confirm
        hr_complexity_filter: if True, reject HR from low-complexity windows

    Returns:
        list of confirmed junction dicts with HR + microhomology classification
    """
    confirmed = []

    for cand in candidates:
        junc_id = cand.get("junction_id", "")
        counts = read_counts.get(junc_id, {})

        spanning = counts.get("spanning_reads", 0)
        total = counts.get("total_reads", 0)
        hq = counts.get("hq_reads", 0)

        if spanning < min_spanning:
            continue

        # Classify junction using single-copy reference analysis
        classification = classify_single_copy_junction(
            ref_seq, cand["start"], cand["end"],
            hr_complexity_filter=hr_complexity_filter,
        )

        # Convert internal 0-based start / Python exclusive end to 1-based inclusive for output
        display_start = cand["start"] + 1
        display_end = cand["end"]

        result = {
            "seq_id": seq_id,
            "dup_start": display_start,
            "dup_end": display_end,
            "dup_size": cand.get("dup_size", cand["end"] - cand["start"]),
            "junction_id": junc_id,
            "junction_sequence": cand["sequence"],
            "junction_k": cand["k"],
            "total_reads": total,
            "hq_reads": hq,
            "spanning_reads": spanning,
            "_screen_spanning": spanning,  # pass 1 count for comparison
            "is_hr_signature": classification["is_hr_signature"],
            "hr_match_len": classification["hr_match_len"],
            "hr_identity": classification["hr_identity"],
            "hr_scenario": classification["hr_scenario"],
            "microhomology_bp": classification["microhomology_bp"],
            "microhomology_seq": classification["microhomology_seq"],
        }

        confirmed.append(result)
        hr_label = "HR" if classification["is_hr_signature"] else f"MH={classification['microhomology_bp']}bp"
        logger.debug(
            f"    Screen hit: {seq_id}:{display_start:,}-{display_end:,} "
            f"({cand['end']-cand['start']:,} bp) "
            f"{hr_label} [{spanning} spanning reads]"
        )

    # Sort by read support
    confirmed.sort(key=lambda x: x["spanning_reads"], reverse=True)

    return confirmed


def _save_results(confirmed, output_dir):
    """Save module 2 results to TSV, JSON, and FASTA."""
    # Remove internal fields
    for junc in confirmed:
        junc.pop("_screen_spanning", None)

    # TSV
    tsv_path = os.path.join(output_dir, "confirmed_junctions.tsv")
    with open(tsv_path, "w") as f:
        if confirmed:
            header = list(confirmed[0].keys())
            f.write("\t".join(header) + "\n")
            for junc in confirmed:
                row = [str(junc.get(h, "")) for h in header]
                f.write("\t".join(row) + "\n")
        else:
            f.write("# No junctions confirmed\n")

    logger.info(f"  Results saved to {tsv_path}")

    # JSON
    json_path = os.path.join(output_dir, "confirmed_junctions.json")
    with open(json_path, "w") as f:
        json.dump(confirmed, f, indent=2)

    # Metadata file for module 3 (if junctions were found)
    if confirmed:
        meta_path = os.path.join(output_dir, "junctions_metadata.tsv")
        with open(meta_path, "w") as f:
            f.write("# Metadata for Tandem module 3 (-pop)\n")
            f.write("# Name\tStart\tEnd\tSequence\n")
            for junc in confirmed:
                hr = junc.get("is_hr_signature", False)
                mh = junc.get("microhomology_bp", 0)
                mech = "HR" if hr else f"MH{mh}bp"
                name = f"{mech}_{junc['dup_start']}_{junc['dup_end']}"
                seq = junc.get("junction_sequence", "")
                f.write(f"{name}\t{junc['dup_start']}\t{junc['dup_end']}\t{seq}\n")
        logger.info(
            f"  Module 3 metadata saved to {meta_path} "
            f"(use with: tandem -r ref.fna -i R1.fq -I R2.fq -pop -m {meta_path})"
        )

        # FASTA of confirmed junction sequences (extended to unique k)
        fasta_path = os.path.join(output_dir, "confirmed_junctions.fasta")
        with open(fasta_path, "w") as f:
            for junc in confirmed:
                hr = junc.get("is_hr_signature", False)
                mh = junc.get("microhomology_bp", 0)
                mech = "HR" if hr else f"MH{mh}bp"
                k = junc.get("junction_k", 0)
                spanning = junc.get("spanning_reads", 0)
                header = (
                    f">{junc['seq_id']}:{junc['dup_start']}-{junc['dup_end']} "
                    f"size={junc.get('dup_size', 0)}bp "
                    f"{mech} k={k} spanning_reads={spanning}"
                )
                f.write(header + "\n")
                seq = junc.get("junction_sequence", "")
                # Write sequence in 80-char lines
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\n")
        logger.info(f"  Confirmed junction FASTA saved to {fasta_path}")
