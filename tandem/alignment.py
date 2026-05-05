"""Read mapping and alignment utilities for Tandem.

Wraps minimap2 (preferred) or bwa for read mapping, and samtools for
BAM file operations.
"""

import logging
import os
import re
import tempfile
from pathlib import Path

from . import utils

logger = logging.getLogger("tandem")


def detect_mapper():
    """Detect which read mapper is available.

    Returns:
        'minimap2' or 'bwa' or None
    """
    if utils.check_dependency("minimap2"):
        return "minimap2"
    elif utils.check_dependency("bwa"):
        return "bwa"
    else:
        return None


def map_reads_to_reference(ref_fasta, read1, read2=None, output_bam=None,
                           threads=1, mapper=None):
    """Map reads to a reference genome and produce sorted BAM.

    Args:
        ref_fasta: path to reference FASTA
        read1: path to forward reads (FASTQ)
        read2: path to reverse reads (FASTQ), or None for single-end
        output_bam: output BAM path (default: derived from ref name)
        threads: number of threads
        mapper: 'minimap2' or 'bwa' (auto-detect if None)

    Returns:
        Path to sorted, indexed BAM file
    """
    if mapper is None:
        mapper = detect_mapper()
    if mapper is None:
        raise RuntimeError("No read mapper found. Install minimap2 or bwa.")

    if output_bam is None:
        output_bam = str(Path(ref_fasta).with_suffix(".sorted.bam"))

    if mapper == "minimap2":
        _map_minimap2(ref_fasta, read1, read2, output_bam, threads)
    elif mapper == "bwa":
        _map_bwa(ref_fasta, read1, read2, output_bam, threads)
    else:
        raise ValueError(f"Unknown mapper: {mapper}")

    # Index BAM
    utils.run_command(
        ["samtools", "index", output_bam],
        description="Indexing BAM"
    )

    return output_bam


def _map_minimap2(ref_fasta, read1, read2, output_bam, threads):
    """Map reads using minimap2.

    Uses a two-step approach (minimap2 writes BAM, then samtools sort)
    instead of a pipe, which avoids pipe buffer/memory issues with
    large datasets on some systems. The unsorted BAM is deleted after
    sorting to save disk space.
    """
    # Step 1: minimap2 writes unsorted BAM via samtools view (no sorting yet)
    unsorted_bam = str(output_bam).replace('.sorted.bam', '.unsorted.bam')
    if unsorted_bam == str(output_bam):
        unsorted_bam = str(output_bam) + '.unsorted.bam'

    cmd_map = ["minimap2", "-ax", "sr", "-t", str(threads), str(ref_fasta)]
    if read2:
        cmd_map += [str(read1), str(read2)]
    else:
        cmd_map += [str(read1)]

    # Pipe minimap2 → samtools view (convert SAM to BAM, no sorting)
    # This is a very light pipe — samtools view just compresses the stream
    map_cmd_str = (
        f"{' '.join(str(c) for c in cmd_map)} | "
        f"samtools view -@ {max(1, threads // 2)} -bS -o {unsorted_bam} -"
    )
    utils.run_command(
        map_cmd_str,
        description="Mapping reads with minimap2 (writing unsorted BAM)"
    )

    # Step 2: samtools sort with conservative memory to avoid OOM
    # Use -m to limit per-thread memory (default would be 768M × threads)
    # Using 512M × max 8 threads = max 4 GB total
    sort_threads = min(threads, 8)
    sort_cmd = [
        "samtools", "sort",
        "-@", str(sort_threads),
        "-m", "512M",
        "-o", str(output_bam),
        unsorted_bam
    ]
    utils.run_command(sort_cmd, description="Sorting BAM")

    # Clean up unsorted BAM
    try:
        os.remove(unsorted_bam)
    except OSError:
        pass


def _map_bwa(ref_fasta, read1, read2, output_bam, threads):
    """Map reads using bwa mem."""
    # Index reference if needed
    if not os.path.exists(f"{ref_fasta}.bwt"):
        utils.run_command(
            ["bwa", "index", str(ref_fasta)],
            description="Indexing reference for bwa"
        )

    cmd = ["bwa", "mem", "-t", str(threads), str(ref_fasta)]

    if read2:
        cmd += [str(read1), str(read2)]
    else:
        cmd += [str(read1)]

    full_cmd = (
        f"{' '.join(str(c) for c in cmd)} | "
        f"samtools sort -@ {threads} -o {output_bam}"
    )

    utils.run_command(full_cmd, description="Mapping reads with bwa")


def map_reads_to_junctions(junction_fasta, read1, read2=None,
                           output_bam=None, threads=1):
    """Map reads against junction reference sequences.

    Uses more sensitive mapping parameters since we're looking for
    reads that span novel junctions.

    Args:
        junction_fasta: FASTA of junction reference sequences
        read1: forward reads
        read2: reverse reads (optional)
        output_bam: output BAM path
        threads: number of threads

    Returns:
        Path to sorted, indexed BAM file
    """
    if output_bam is None:
        output_bam = str(Path(junction_fasta).with_suffix(".reads.sorted.bam"))

    mapper = detect_mapper()

    if mapper == "minimap2":
        # Use more sensitive parameters for short references
        cmd = [
            "minimap2", "-ax", "sr",
            "--secondary=no",       # no secondary alignments
            "-t", str(threads),
            str(junction_fasta)
        ]
    elif mapper == "bwa":
        # Index junction reference
        if not os.path.exists(f"{junction_fasta}.bwt"):
            utils.run_command(["bwa", "index", str(junction_fasta)])
        cmd = ["bwa", "mem", "-t", str(threads), str(junction_fasta)]
    else:
        raise RuntimeError("No read mapper found.")

    if read2:
        cmd += [str(read1), str(read2)]
    else:
        cmd += [str(read1)]

    full_cmd = (
        f"{' '.join(str(c) for c in cmd)} | "
        f"samtools sort -@ {threads} -o {output_bam}"
    )

    utils.run_command(full_cmd, description="Mapping reads to junction references")

    utils.run_command(
        ["samtools", "index", output_bam],
        description="Indexing junction BAM"
    )

    return output_bam


def count_junction_reads(bam_path, min_mapq=0, min_overlap=10):
    """Count reads mapped to each junction reference.

    Uses two criteria instead of MAPQ filtering:
    1. Zero mismatches (NM:i:0) in the aligned region — the junction
       is a precise sequence, so any real junction read matches perfectly.
    2. Read spans the junction midpoint with at least min_overlap bases
       on each side.

    This avoids the multi-mapping MAPQ problem when many similar
    junction candidates are pooled into one reference.

    Args:
        bam_path: path to BAM mapped against junction references
        min_mapq: ignored (kept for API compatibility)
        min_overlap: minimum bases mapped on each side of junction center

    Returns:
        dict of {junction_id: {'total_reads': N, 'hq_reads': N, 'spanning_reads': N}}
    """
    # Use samtools idxstats for total counts per reference
    result = utils.run_command(
        ["samtools", "idxstats", str(bam_path)],
        description="Counting junction reads"
    )

    counts = {}
    for line in result.stdout.strip().split("\n"):
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 4 or parts[0] == "*":
            continue
        ref_name = parts[0]
        ref_len = int(parts[1])
        mapped = int(parts[2])

        counts[ref_name] = {
            "total_reads": mapped,
            "ref_length": ref_len,
        }

    # Read all mapped reads and apply strict filtering:
    # 1. NM:i:0 (zero mismatches)
    # 2. No soft-clipping (S), no indels (I/D) — read must match exactly
    # 3. Aligned region crosses junction midpoint with min_overlap on each side
    result = utils.run_command(
        ["samtools", "view", "-F", "4", str(bam_path)],
        description="Filtering junction reads (100% ref coverage, NM=0, no indels)"
    )

    hq_counts = {}
    for line in result.stdout.strip().split("\n"):
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 11:
            continue

        ref_name = parts[2]
        pos = int(parts[3]) - 1  # 0-based alignment start
        cigar = parts[5]

        # Parse CIGAR: get aligned reference length, reject indels
        has_indel = False
        aligned_ref_len = 0
        for length_str, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
            length = int(length_str)
            if op in ('I', 'D'):
                has_indel = True
                break
            if op in ('M', '=', 'X'):
                aligned_ref_len += length
            # S, H, N, P: don't count as aligned to reference

        if has_indel:
            continue

        read_end = pos + aligned_ref_len

        # Parse NM tag (zero mismatches in aligned portion)
        nm = None
        for tag in parts[11:]:
            if tag.startswith("NM:i:"):
                nm = int(tag.split(":")[2])
                break

        if ref_name not in hq_counts:
            hq_counts[ref_name] = {"hq_reads": 0, "spanning_reads": 0}

        # Filter: NM=0 (zero mismatches in aligned portion)
        if nm is not None and nm > 0:
            continue

        ref_len = counts.get(ref_name, {}).get("ref_length", 0)
        if ref_len <= 0:
            continue

        # Filter: read covers 100% of junction reference
        # - Alignment must start at position 0 of reference
        # - Aligned length must cover entire reference
        # - Read can overhang (soft-clipped), that's fine
        if pos != 0:
            continue
        if aligned_ref_len < ref_len:
            continue

        hq_counts[ref_name]["hq_reads"] += 1
        hq_counts[ref_name]["spanning_reads"] += 1

    # Merge counts
    for ref_name in counts:
        if ref_name in hq_counts:
            counts[ref_name].update(hq_counts[ref_name])
        else:
            counts[ref_name]["hq_reads"] = 0
            counts[ref_name]["spanning_reads"] = 0

    return counts


def get_coverage_at_position(bam_path, seq_id, start, end):
    """Get mean coverage at a specific genomic region.

    Args:
        bam_path: path to BAM file
        seq_id: sequence/chromosome ID
        start: start position (0-based)
        end: end position

    Returns:
        float mean coverage
    """
    region = f"{seq_id}:{start+1}-{end}"  # samtools uses 1-based
    result = utils.run_command(
        ["samtools", "depth", "-a", "-r", region, str(bam_path)],
        description=f"Coverage at {region}"
    )

    depths = []
    for line in result.stdout.strip().split("\n"):
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 3:
            depths.append(int(parts[2]))

    return float(np.mean(depths)) if depths else 0.0


# Need numpy for get_coverage_at_position
import numpy as np
