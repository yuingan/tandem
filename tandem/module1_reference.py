"""Module 1: Reference genome tandem duplication detection.

Three-stage pipeline for detecting and characterizing tandem duplications
in a reference genome:

    Stage 1 — Detection (sequence-geometric):
        NUCmer self-alignment finds candidate tandem duplications.
        Pure sequence-based; does not depend on annotation.

    Stage 2 — CDS subcategorization (optional, annotation-based):
        If a GFF or GenBank annotation is provided, each tandem is tagged
        with CDS content statistics (n_cds_copy1, n_cds_gap, etc.) and a
        content_category (intergenic / tandem_single_gene / tandem_segmental /
        proximal_*). Does not affect detection or mechanism inference.

    Stage 3 — Mechanism classification (sequence-only):
        HR signature detection at outer boundaries + microhomology bp at the
        inner junction. Reports mechanism_confidence (low/moderate/high)
        based on copy length, since junction-based inference is less reliable
        at short lengths.

Supports circular genomes: distance between copies is computed as the
minimum of linear distance and wrap-around distance through the origin.

Coordinate system:
    Coordinates throughout are 1-based inclusive (paper/NCBI/GenBank/IGV
    convention). NUCmer show-coords outputs 1-based inclusive, so we store
    them as-is. For Python slicing (classification), we convert internally:
        0-based start = 1-based start - 1
        0-based end   = 1-based end  (inclusive -> exclusive)

Usage: tandem -r reference.fna
       tandem -r reference.fna --annotation annotation.gff
"""

import logging
import os
import json
from pathlib import Path

from . import utils
from .junction import KmerIndex, classify_junction, HR_MIN_CONSECUTIVE

logger = logging.getLogger("tandem")


# Mechanism confidence thresholds (copy length in bp)
CONFIDENCE_LOW_MAX = 300      # <300bp: flanking windows overlap copy itself
CONFIDENCE_MODERATE_MAX = 500 # 300-500bp: moderate
# >=500bp: high


def run_nucmer_self_alignment(ref_fasta, output_dir, threads=1,
                               min_match=20, min_cluster=65):
    """Run NUCmer self-alignment pipeline.

    Returns: Path to coords file
    """
    prefix = os.path.join(output_dir, "self_align")
    delta_file = prefix + ".delta"
    filtered_delta = prefix + ".filtered.delta"
    coords_file = prefix + ".coords"

    cmd = [
        "nucmer",
        "--maxmatch",
        "-l", str(min_match),
        "-c", str(min_cluster),
        "-t", str(threads),
        "--prefix", prefix,
        str(ref_fasta),
        str(ref_fasta),
    ]
    utils.run_command(cmd, description="Running NUCmer self-alignment")

    cmd = f"delta-filter -i 80 -l 50 {delta_file}"
    result = utils.run_command(cmd, description="Filtering alignments")
    with open(filtered_delta, "w") as f:
        f.write(result.stdout)

    cmd = f"show-coords -r -c -l -T {filtered_delta}"
    result = utils.run_command(cmd, description="Extracting coordinates")
    with open(coords_file, "w") as f:
        f.write(result.stdout)

    return coords_file


def parse_coords(coords_file):
    """Parse show-coords tab-delimited output.

    show-coords -r -c -l -T produces 13 columns:
    [0] S1  [1] E1  [2] S2  [3] E2  (1-based inclusive)
    [4] LEN1  [5] LEN2  [6] %IDY
    [7] LEN_R  [8] LEN_Q  [9] COV_R  [10] COV_Q
    [11] REF_TAG  [12] QRY_TAG
    """
    alignments = []

    with open(coords_file) as f:
        lines = f.readlines()

    start_idx = 0
    for i, line in enumerate(lines):
        if '[S1]' in line:
            start_idx = i + 1
            break

    for line in lines[start_idx:]:
        parts = line.strip().split('\t')
        if len(parts) != 13:
            continue

        try:
            alignments.append({
                'ref_start': int(parts[0]),
                'ref_end': int(parts[1]),
                'qry_start': int(parts[2]),
                'qry_end': int(parts[3]),
                'ref_alen': int(parts[4]),
                'qry_alen': int(parts[5]),
                'identity': float(parts[6]),
                'ref_id': parts[11],
                'qry_id': parts[12],
            })
        except (ValueError, IndexError):
            continue

    return alignments


SHORT_DUP_THRESHOLD = 500        # bp - below this, apply extra filters
COMPLEXITY_MIN_FRACTION = 0.20   # <500bp: at least 20% of possible 4-mers present
COMPLEXITY_KMER = 4              # size of k-mer for complexity calculation


def _sequence_complexity_ok(seq, k=COMPLEXITY_KMER,
                             min_fraction=COMPLEXITY_MIN_FRACTION):
    """Quick complexity check: fraction of possible k-mers present in seq.

    A mono/dinucleotide-repeat sequence (AAAAAA... or ATATATAT...) will
    contain only 1 unique k-mer, while random sequence will contain many.
    This is a poor-man's low-complexity filter that avoids needing
    dustmasker or another external dependency.

    Returns True if the sequence has adequate complexity, False if it
    is dominated by repeats of a few bases.
    """
    if len(seq) < k:
        # Too short to evaluate — accept (don't filter on length)
        return True

    kmers = set()
    for i in range(len(seq) - k + 1):
        kmers.add(seq[i : i + k])

    possible = 4 ** k  # 256 for k=4
    # Normalize by min(possible, n_kmers_total) so short sequences aren't penalized
    max_possible = min(possible, len(seq) - k + 1)
    if max_possible <= 0:
        return True
    return (len(kmers) / max_possible) >= min_fraction


def identify_tandem_duplications(alignments, seq_lengths, max_distance=50000,
                                  min_identity=80.0, min_size=200,
                                  circular=True, sequences=None,
                                  short_dup_threshold=SHORT_DUP_THRESHOLD,
                                  apply_short_dup_filters=True):
    """Identify tandem duplications from NUCmer self-alignment results.

    Short-duplication filters (applied to duplications below
    short_dup_threshold, default 500bp, when apply_short_dup_filters=True):

      1. Gap < copy length: the inter-copy gap must be shorter than the
         copy itself. Prevents calling chance short matches that happen
         to be slightly closer than max_distance.

      2. Sequence complexity: the duplicated copy must have at least
         COMPLEXITY_MIN_FRACTION of possible 4-mers present, excluding
         mono/dinucleotide-repeat matches (e.g. AAAAAA, ATATATAT).

    These filters reduce false positives for short duplications, where
    random sequence similarity is more likely to produce a signal.

    Args:
        alignments: list of alignment dicts from parse_coords
        seq_lengths: dict of {seq_id: length} for circular distance
        max_distance: maximum distance between copies to consider tandem
        min_identity: minimum alignment identity
        min_size: minimum duplication size (default 200)
        circular: if True, compute circular distance
        sequences: dict of {seq_id: sequence} — required if complexity
                   check is applied (i.e. when short-dup filters are on)
        short_dup_threshold: bp threshold below which extra filters apply
        apply_short_dup_filters: if False, skip short-duplication filters

    Returns:
        list of tandem duplication dicts
    """
    tandems = []
    seen = set()
    n_filtered_by_short = 0
    n_filtered_by_complexity = 0

    for aln in alignments:
        # Same contig only
        if aln["ref_id"] != aln["qry_id"]:
            continue

        if aln["identity"] < min_identity:
            continue
        if aln["ref_alen"] < min_size:
            continue

        # Normalize positions (1-based inclusive)
        pos1 = sorted([aln["ref_start"], aln["ref_end"]])
        pos2 = sorted([aln["qry_start"], aln["qry_end"]])

        # Skip self-match (nearly identical coordinates)
        if abs(pos1[0] - pos2[0]) < 10 and abs(pos1[1] - pos2[1]) < 10:
            continue

        # Calculate distance and assign copy order
        nucmer_overlap_bp = 0

        if pos1[1] < pos2[0]:
            linear_distance = pos2[0] - pos1[1] - 1
            copy1, copy2 = pos1, pos2
        elif pos2[1] < pos1[0]:
            linear_distance = pos1[0] - pos2[1] - 1
            copy1, copy2 = pos2, pos1
        else:
            # Overlapping copies
            overlap = min(pos1[1], pos2[1]) - max(pos1[0], pos2[0]) + 1
            if overlap > aln["ref_alen"] * 0.1:
                continue

            nucmer_overlap_bp = overlap
            if pos1[0] <= pos2[0]:
                copy1, copy2 = list(pos1), list(pos2)
            else:
                copy1, copy2 = list(pos2), list(pos1)

            copy2[0] = copy1[1] + 1
            linear_distance = 0

        # Circular distance: for circular genomes, the true distance
        # may be shorter going around the origin
        seq_id = aln["ref_id"]
        genome_len = seq_lengths.get(seq_id, 0)

        if circular and genome_len > 0 and linear_distance > 0:
            # Wrap-around distance through the origin
            # For two copies at positions copy1 and copy2 on a circular genome:
            # circular_distance = genome_length - (copy2_end - copy1_start)
            span = copy2[1] - copy1[0]
            wrap_distance = genome_len - span
            distance = min(linear_distance, max(0, wrap_distance))
            is_circular_closer = wrap_distance < linear_distance
        else:
            distance = linear_distance
            is_circular_closer = False

        if distance > max_distance:
            continue

        # Determine orientation
        ref_fwd = aln["ref_start"] < aln["ref_end"]
        qry_fwd = aln["qry_start"] < aln["qry_end"]
        orientation = "direct" if ref_fwd == qry_fwd else "inverted"

        # Classify proximity
        if nucmer_overlap_bp > 0:
            proximity = "adjacent"
        elif distance <= 0:
            proximity = "adjacent"
        elif distance <= 500:
            proximity = "adjacent"
        elif distance <= 2000:
            proximity = "near"
        else:
            proximity = "proximal"

        # Deduplicate
        key = (seq_id, copy1[0], copy1[1], copy2[0], copy2[1])
        if key in seen:
            continue
        seen.add(key)

        # Short-duplication filters (<short_dup_threshold bp)
        dup_size = aln["ref_alen"]
        if apply_short_dup_filters and dup_size < short_dup_threshold:
            # Filter 1: gap must be shorter than copy length
            if distance >= dup_size:
                n_filtered_by_short += 1
                continue

            # Filter 2: complexity check (if sequences provided)
            if sequences is not None:
                ref = sequences.get(seq_id, "")
                if ref:
                    # Extract copy 1 sequence (1-based inclusive coords)
                    copy_seq = ref[copy1[0] - 1 : copy1[1]]
                    if not _sequence_complexity_ok(copy_seq):
                        n_filtered_by_complexity += 1
                        continue

        tandems.append({
            "seq_id": seq_id,
            "copy1_start": copy1[0],
            "copy1_end": copy1[1],
            "copy2_start": copy2[0],
            "copy2_end": copy2[1],
            "distance": distance,
            "nucmer_overlap_bp": nucmer_overlap_bp,
            "circular_closer": is_circular_closer,
            "size": aln["ref_alen"],
            "identity": aln["identity"],
            "orientation": orientation,
            "proximity_class": proximity,
            "distance_to_size_ratio": round(distance / aln["ref_alen"], 2) if aln["ref_alen"] > 0 else 0.0,
            "detection_source": "native",
        })

    tandems.sort(key=lambda t: (t["seq_id"], t["copy1_start"]))

    # Count circular-closer hits
    n_circular = sum(1 for t in tandems if t["circular_closer"])
    if n_circular > 0:
        logger.info(f"  {n_circular} tandem(s) closer via circular distance (through origin)")

    # Log short-dup filter stats
    if n_filtered_by_short > 0 or n_filtered_by_complexity > 0:
        logger.info(
            f"  Short-duplication filters removed: "
            f"{n_filtered_by_short} for gap >= copy length, "
            f"{n_filtered_by_complexity} for low complexity"
        )

    logger.info(f"  Identified {len(tandems)} tandem duplication candidates")
    return tandems


def classify_tandem_mechanisms(tandems, sequences, flag_distance=2000,
                                hr_min_consec=None, hr_min_identity=None,
                                hr_max_inward=None, hr_outward_ext=None,
                                hr_complexity_filter=True):
    """Classify the mechanism for each tandem duplication.

    Uses the two-step approach:
    1. HR signature detection (flanking repeat alignment)
    2. Microhomology bp measurement

    Classification notes (flags):
    - 'junction_distant': distance > flag_distance, junction signal may be degraded
    - 'too_small_for_HR_check': copy too small for meaningful HR detection

    Args:
        tandems: list of tandem duplication dicts (1-based coords)
        sequences: dict of {seq_id: sequence_string}
        flag_distance: distance threshold for junction_distant flag
        hr_min_consec: override for HR_MIN_CONSECUTIVE (default: module default)
        hr_min_identity: override for HR_MIN_IDENTITY (default: module default)
        hr_max_inward: override for HR_MAX_WINDOW (default: module default)
        hr_outward_ext: override for HR_OUTWARD_EXTENSION (default: module default)
        hr_complexity_filter: if True, reject HR calls where both flanking
            windows are low-complexity (microsatellite periodicity). Default True.

    Returns:
        Updated list with classification fields added
    """
    # Build HR kwargs, falling back to module defaults for any not provided
    from .junction import (
        HR_MIN_CONSECUTIVE as DEF_HR_MIN_CONSEC,
        HR_MIN_IDENTITY as DEF_HR_MIN_ID,
        HR_MAX_WINDOW as DEF_HR_MAX_WIN,
        HR_OUTWARD_EXTENSION as DEF_HR_OUT_EXT,
    )
    hr_kwargs = {
        "hr_min_consec": hr_min_consec if hr_min_consec is not None else DEF_HR_MIN_CONSEC,
        "hr_min_identity": hr_min_identity if hr_min_identity is not None else DEF_HR_MIN_ID,
        "hr_max_inward": hr_max_inward if hr_max_inward is not None else DEF_HR_MAX_WIN,
        "hr_outward_ext": hr_outward_ext if hr_outward_ext is not None else DEF_HR_OUT_EXT,
        "hr_complexity_filter": hr_complexity_filter,
    }

    filter_label = "ON" if hr_complexity_filter else "OFF"
    logger.info(
        f"  HR detection params: min_consec={hr_kwargs['hr_min_consec']}bp, "
        f"min_identity={hr_kwargs['hr_min_identity']}, "
        f"max_inward={hr_kwargs['hr_max_inward']}bp, "
        f"outward_ext={hr_kwargs['hr_outward_ext']}bp, "
        f"complexity_filter={filter_label}"
    )

    classified = []
    hr_count = 0
    total = 0

    for td in tandems:
        seq = sequences.get(td["seq_id"])
        if seq is None:
            logger.warning(f"  Sequence {td['seq_id']} not found, skipping")
            continue

        # Convert 1-based inclusive -> 0-based exclusive
        c1s_0 = td["copy1_start"] - 1
        c1e_0 = td["copy1_end"]
        c2s_0 = td["copy2_start"] - 1
        c2e_0 = td["copy2_end"]

        # Classify junction
        classification = classify_junction(seq, c1s_0, c1e_0, c2s_0, c2e_0, **hr_kwargs)

        td["is_hr_signature"] = classification["is_hr_signature"]
        td["hr_match_len"] = classification["hr_match_len"]
        td["hr_identity"] = classification["hr_identity"]
        td["hr_scenario"] = classification["hr_scenario"]
        td["microhomology_bp"] = classification["microhomology_bp"]
        td["microhomology_seq"] = classification["microhomology_seq"]
        td["junction_gap_bp"] = classification["junction_gap_bp"]
        td["junction_seq"] = classification["junction_seq"]

        # NUCmer overlap → microhomology consolidation.
        # When NUCmer reports overlapping copies (nucmer_overlap_bp > 0),
        # the overlap resolution shifts copy2's start past copy1's end.
        # find_microhomology then sees abutting copies and returns 0.
        # But the NUCmer overlap IS the microhomology — it's the shared
        # sequence at the junction that NUCmer detected as alignment overlap.
        # Consolidate: use the larger of the two values.
        nucmer_overlap = td.get("nucmer_overlap_bp", 0)
        if nucmer_overlap > td["microhomology_bp"]:
            td["microhomology_bp"] = nucmer_overlap
            # Can't recover exact sequence after overlap resolution
            td["microhomology_seq"] = f"(from_nucmer_overlap_{nucmer_overlap}bp)"

        # Mechanism confidence based on copy length.
        # Below 300bp, the adaptive HR search window (dup_length // 2)
        # is small relative to the copy, making mechanism inference
        # less reliable. Above 500bp, flanking inference is robust.
        copy_len = td["size"]
        if copy_len < CONFIDENCE_LOW_MAX:
            td["mechanism_confidence"] = "low"
        elif copy_len < CONFIDENCE_MODERATE_MAX:
            td["mechanism_confidence"] = "moderate"
        else:
            td["mechanism_confidence"] = "high"

        # Build classification_note flags
        notes = []
        if td["distance"] > flag_distance:
            notes.append("junction_distant")
        if td["size"] < 2 * HR_MIN_CONSECUTIVE:
            notes.append("too_small_for_HR_check")

        td["classification_note"] = ";".join(notes) if notes else ""

        classified.append(td)
        total += 1
        if classification["is_hr_signature"]:
            hr_count += 1

    logger.info(
        f"  Classification: {total} total, "
        f"{hr_count} HR-like, {total - hr_count} non-HR"
    )

    # Log flag counts
    flag_counts = {}
    for td in classified:
        for note in td["classification_note"].split(";"):
            if note:
                flag_counts[note] = flag_counts.get(note, 0) + 1
    if flag_counts:
        logger.info(f"  Classification flags:")
        for flag, count in sorted(flag_counts.items(), key=lambda x: -x[1]):
            logger.info(f"    {flag}: {count}")

    # Log microhomology distribution for non-HR
    mh_dist = {}
    for td in classified:
        if not td["is_hr_signature"]:
            mh = td["microhomology_bp"]
            mh_dist[mh] = mh_dist.get(mh, 0) + 1
    if mh_dist:
        top_mh = sorted(mh_dist.items(), key=lambda x: -x[1])[:10]
        logger.info(f"  Non-HR microhomology distribution (top 10):")
        for mh_bp, count in top_mh:
            logger.info(f"    {mh_bp} bp: {count}")

    return classified


def _parse_gff_cds_features(annotation_path):
    """Parse a GFF3 file and return a list of CDS features.

    Returns:
        dict of seq_id -> list of (start, end, strand) tuples in 1-based
        inclusive coordinates.
    """
    cds_by_seq = {}
    with open(annotation_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 8:
                continue
            feature_type = parts[2].lower()
            if feature_type != 'cds':
                continue
            seq_id = parts[0]
            try:
                start = int(parts[3])  # 1-based inclusive
                end = int(parts[4])
            except ValueError:
                continue
            strand = parts[6]
            if seq_id not in cds_by_seq:
                cds_by_seq[seq_id] = []
            cds_by_seq[seq_id].append((start, end, strand))

    # Sort by start position for each sequence
    for seq_id in cds_by_seq:
        cds_by_seq[seq_id].sort(key=lambda x: x[0])

    return cds_by_seq


def _parse_genbank_cds_features(annotation_path):
    """Parse a GenBank file and return a list of CDS features.

    Returns:
        dict of seq_id -> list of (start, end, strand) tuples in 1-based
        inclusive coordinates.
    """
    try:
        from Bio import SeqIO
    except ImportError:
        logger.error(
            "Biopython required to parse GenBank annotations. "
            "Install with: pip install biopython"
        )
        return {}

    cds_by_seq = {}
    for record in SeqIO.parse(annotation_path, "genbank"):
        seq_id = record.id
        cds_by_seq[seq_id] = []
        for feature in record.features:
            if feature.type.lower() != 'cds':
                continue
            # Biopython uses 0-based half-open; convert to 1-based inclusive
            start_1based = int(feature.location.start) + 1
            end_1based = int(feature.location.end)
            strand = '+' if feature.location.strand == 1 else '-'
            cds_by_seq[seq_id].append((start_1based, end_1based, strand))

    for seq_id in cds_by_seq:
        cds_by_seq[seq_id].sort(key=lambda x: x[0])

    return cds_by_seq


def _count_cds_in_interval(cds_list, interval_start, interval_end,
                            min_overlap_frac=0.5):
    """Count CDS fully contained (or >= min_overlap_frac overlapping) in interval.

    Args:
        cds_list: list of (cds_start, cds_end, strand) tuples (1-based inclusive)
        interval_start: 1-based inclusive start
        interval_end: 1-based inclusive end
        min_overlap_frac: fraction of CDS that must be inside the interval

    Returns:
        (n_fully_contained, n_partial)
    """
    n_full = 0
    n_partial = 0
    for cds_start, cds_end, _ in cds_list:
        # Skip CDSs entirely outside the interval
        if cds_end < interval_start or cds_start > interval_end:
            continue

        cds_length = cds_end - cds_start + 1
        if cds_length <= 0:
            continue

        overlap_start = max(cds_start, interval_start)
        overlap_end = min(cds_end, interval_end)
        overlap = overlap_end - overlap_start + 1
        overlap_frac = overlap / cds_length

        if cds_start >= interval_start and cds_end <= interval_end:
            n_full += 1
        elif overlap_frac >= min_overlap_frac:
            n_partial += 1

    return n_full, n_partial


def _classify_content_category(n_cds_copy1, n_cds_gap, gap_bp):
    """Classify duplication by its CDS content.

    Categories:
        intergenic — no CDSs in either copy
        tandem_single_gene — 1 CDS per copy, 0 intervening CDSs
        tandem_segmental — >1 CDS per copy, 0 intervening CDSs
        proximal_single_gene — 1 CDS per copy, 1-2 intervening CDSs
        proximal_segmental — >1 CDS per copy, 1-2 intervening CDSs
        other — larger gap or mixed composition
    """
    if n_cds_copy1 == 0:
        return "intergenic"

    if n_cds_gap == 0:
        # True tandem (copies abutting)
        if n_cds_copy1 == 1:
            return "tandem_single_gene"
        else:
            return "tandem_segmental"
    elif n_cds_gap <= 2:
        # Proximal (small gap with a few intervening CDSs)
        if n_cds_copy1 == 1:
            return "proximal_single_gene"
        else:
            return "proximal_segmental"
    else:
        return "other"


def _classify_boundary_type(cds_list, boundary_pos):
    """Determine whether a boundary falls inside a CDS, at a CDS edge, or in IGS.

    Args:
        cds_list: list of (cds_start, cds_end, strand) tuples (1-based inclusive)
        boundary_pos: 1-based coordinate to check

    Returns:
        str: "intragenic" | "cds_boundary" | "intergenic"
    """
    EDGE_TOLERANCE = 3  # bp tolerance for "at CDS edge"
    for cds_start, cds_end, _ in cds_list:
        if boundary_pos < cds_start - EDGE_TOLERANCE:
            # We've passed the boundary (sorted list)
            return "intergenic"
        # Check if at start edge
        if abs(boundary_pos - cds_start) <= EDGE_TOLERANCE:
            return "cds_boundary"
        # Check if at end edge
        if abs(boundary_pos - cds_end) <= EDGE_TOLERANCE:
            return "cds_boundary"
        # Check if inside CDS
        if cds_start <= boundary_pos <= cds_end:
            return "intragenic"
    return "intergenic"


def subcategorize_by_annotation(tandems, annotation_path):
    """Stage 2: annotate tandems with CDS content categories.

    For each tandem duplication, count CDSs fully contained in copy 1,
    copy 2, and the intervening gap region. Assign a content_category
    (intergenic / tandem_single_gene / tandem_segmental / proximal_*) and
    a boundary_type describing where the duplication boundaries sit
    relative to CDS features.

    Args:
        tandems: list of tandem dicts from Stage 1 (with 1-based copy1_start,
                 copy1_end, copy2_start, copy2_end)
        annotation_path: path to GFF3 or GenBank file

    Returns:
        Updated list with fields added: n_cds_copy1, n_cds_copy2, n_cds_gap,
        content_category, boundary_type. If annotation parsing fails, each
        tandem gets "NA" for all Stage 2 fields.
    """
    # Detect format from extension
    path = str(annotation_path).lower()
    if path.endswith('.gff') or path.endswith('.gff3'):
        cds_by_seq = _parse_gff_cds_features(annotation_path)
    elif path.endswith('.gbk') or path.endswith('.gb') or path.endswith('.genbank'):
        cds_by_seq = _parse_genbank_cds_features(annotation_path)
    else:
        logger.warning(
            f"  Cannot infer annotation format from '{annotation_path}'. "
            f"Expected .gff, .gff3, .gbk, .gb, or .genbank. "
            f"Skipping Stage 2 subcategorization."
        )
        for td in tandems:
            td["n_cds_copy1"] = "NA"
            td["n_cds_copy2"] = "NA"
            td["n_cds_gap"] = "NA"
            td["content_category"] = "NA"
            td["boundary_type"] = "NA"
        return tandems

    total_cds = sum(len(v) for v in cds_by_seq.values())
    logger.info(
        f"  Parsed {total_cds:,} CDS features across "
        f"{len(cds_by_seq)} sequence(s) from annotation"
    )

    category_counts = {}

    for td in tandems:
        seq_id = td["seq_id"]
        cds_list = cds_by_seq.get(seq_id, [])

        c1_start = td["copy1_start"]
        c1_end = td["copy1_end"]
        c2_start = td["copy2_start"]
        c2_end = td["copy2_end"]

        # Count CDSs in each copy
        n_c1_full, _ = _count_cds_in_interval(cds_list, c1_start, c1_end)
        n_c2_full, _ = _count_cds_in_interval(cds_list, c2_start, c2_end)

        # Count CDSs in gap (between copies)
        gap_bp = td.get("distance", 0)
        if gap_bp > 0 and c2_start > c1_end:
            n_gap_full, _ = _count_cds_in_interval(
                cds_list, c1_end + 1, c2_start - 1
            )
        else:
            n_gap_full = 0

        # Classify content category
        content_category = _classify_content_category(
            n_c1_full, n_gap_full, gap_bp
        )

        # Classify boundary type (use copy 1 start and copy 2 end as outer boundaries)
        b_start = _classify_boundary_type(cds_list, c1_start)
        b_end = _classify_boundary_type(cds_list, c2_end)
        if b_start == b_end:
            boundary_type = b_start
        else:
            boundary_type = f"{b_start}/{b_end}"

        td["n_cds_copy1"] = n_c1_full
        td["n_cds_copy2"] = n_c2_full
        td["n_cds_gap"] = n_gap_full
        td["content_category"] = content_category
        td["boundary_type"] = boundary_type

        category_counts[content_category] = category_counts.get(
            content_category, 0) + 1

    # Log summary
    logger.info("  Content categories:")
    for cat, count in sorted(category_counts.items(),
                             key=lambda x: -x[1]):
        logger.info(f"    {cat}: {count}")

    return tandems


def _load_detection_input_tsv(tsv_path, sequences, seq_lengths,
                               circular=True, flag_distance=2000):
    """Load pre-computed tandem duplication coordinates from a TSV file.

    This lets users feed tandem coordinates from an external tool
    (e.g. SegMantX, BISER, breseq, manual curation) and skip the native
    NUCmer detection. Tandem then runs Stages 2 and 3 on the provided
    coordinates.

    Expected TSV format (tab-separated):
        seq_id  copy1_start  copy1_end  copy2_start  copy2_end
            [size] [identity] [orientation]

    Coordinates are 1-based inclusive (matching tandem's convention).
    Lines starting with '#' are treated as comments. The first
    non-comment row may be a header (auto-detected if its copy1_start
    is not numeric).

    Returns:
        list of tandem duplication dicts compatible with Stage 3
    """
    tandems = []
    n_skipped = 0

    with open(tsv_path) as f:
        header_skipped = False
        for line_num, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                n_skipped += 1
                logger.warning(
                    f"  Line {line_num}: expected >=5 columns, got {len(parts)}. "
                    f"Skipping."
                )
                continue

            # Try to detect header row (first row with non-numeric start)
            if not header_skipped:
                try:
                    int(parts[1])
                except ValueError:
                    logger.info(f"  Treating line {line_num} as header, skipping")
                    header_skipped = True
                    continue
                header_skipped = True

            try:
                seq_id = parts[0]
                copy1_start = int(parts[1])
                copy1_end = int(parts[2])
                copy2_start = int(parts[3])
                copy2_end = int(parts[4])
            except ValueError:
                n_skipped += 1
                continue

            # Optional columns
            size_val = None
            identity_val = 100.0
            orientation = "direct"
            if len(parts) > 5:
                try:
                    size_val = int(parts[5])
                except (ValueError, IndexError):
                    pass
            if len(parts) > 6:
                try:
                    identity_val = float(parts[6])
                except (ValueError, IndexError):
                    pass
            if len(parts) > 7:
                orientation = parts[7].strip().lower() or "direct"

            # Validate seq_id is in reference
            if seq_id not in seq_lengths:
                logger.warning(
                    f"  Line {line_num}: seq_id '{seq_id}' not in reference. "
                    f"Skipping."
                )
                n_skipped += 1
                continue

            if size_val is None:
                size_val = copy1_end - copy1_start + 1

            # Distance (gap) between copies
            if copy2_start > copy1_end:
                distance = copy2_start - copy1_end - 1
            else:
                distance = 0

            # Circular distance
            genome_len = seq_lengths[seq_id]
            if circular and genome_len > 0 and distance > 0:
                span = copy2_end - copy1_start
                wrap_distance = genome_len - span
                is_circular_closer = wrap_distance < distance
                distance = min(distance, max(0, wrap_distance))
            else:
                is_circular_closer = False

            # Proximity class
            if distance <= 500:
                proximity = "adjacent"
            elif distance <= 2000:
                proximity = "near"
            else:
                proximity = "proximal"

            tandems.append({
                "seq_id": seq_id,
                "copy1_start": copy1_start,
                "copy1_end": copy1_end,
                "copy2_start": copy2_start,
                "copy2_end": copy2_end,
                "distance": distance,
                "nucmer_overlap_bp": 0,
                "circular_closer": is_circular_closer,
                "size": size_val,
                "identity": identity_val,
                "orientation": orientation,
                "proximity_class": proximity,
                "distance_to_size_ratio": round(distance / size_val, 2) if size_val > 0 else 0.0,
                "detection_source": "external",
            })

    logger.info(
        f"  Loaded {len(tandems)} tandem duplications from {tsv_path}"
    )
    if n_skipped > 0:
        logger.info(f"  ({n_skipped} lines skipped)")

    return tandems


def run_module1(ref_fasta, output_dir, threads=1, max_distance=50000,
                min_identity=80.0, min_size=200, circular=True,
                flag_distance=2000, annotation=None,
                apply_short_dup_filters=True,
                detection_input=None,
                hr_min_consec=None, hr_min_identity=None,
                hr_max_inward=None, hr_outward_ext=None,
                hr_complexity_filter=True):
    """Run the complete three-stage module 1 pipeline.

    Stage 1 — Detection (sequence-geometric, always runs):
      NUCmer self-alignment + tandem filtering.
      Alternatively, coordinates can be loaded from an external TSV
      (e.g. SegMantX output) via the `detection_input` argument,
      bypassing NUCmer entirely.

    Stage 2 — CDS subcategorization (runs only if annotation is provided):
      Count CDSs in each copy and in the gap region, assign content category.

    Stage 3 — Mechanism classification (sequence-only, always runs):
      HR signature detection + microhomology bp + mechanism confidence.

    Args:
        ref_fasta: path to reference FASTA
        output_dir: output directory
        threads: number of threads
        max_distance: maximum distance between copies (default 50000)
        min_identity: minimum alignment identity (default 80.0)
        min_size: minimum duplication size in bp (default 200)
        circular: compute circular distance for bacterial genomes
        flag_distance: distance threshold for junction_distant flag
        annotation: optional path to GFF3 or GenBank file (enables Stage 2)
        apply_short_dup_filters: apply extra filters for duplications
            below SHORT_DUP_THRESHOLD bp (default True)
        detection_input: optional TSV with pre-computed coordinates
            (columns: seq_id, copy1_start, copy1_end, copy2_start, copy2_end;
            optional: size, identity, orientation). If provided, Stage 1
            NUCmer detection is skipped.
        hr_min_consec: override HR_MIN_CONSECUTIVE
        hr_min_identity: override HR_MIN_IDENTITY
        hr_max_inward: override HR_MAX_WINDOW
        hr_outward_ext: override HR_OUTWARD_EXTENSION

    Returns:
        list of classified tandem duplications
    """
    logger.info("=" * 60)
    logger.info("Module 1: Reference genome tandem duplication detection")
    logger.info("=" * 60)

    output_dir = utils.ensure_dir(output_dir)

    # Validate file extension
    utils.validate_file_extension(ref_fasta, utils.REFERENCE_EXTENSIONS, "Reference")

    # Load reference (handles .gz natively)
    logger.info(f"Loading reference: {ref_fasta}")
    sequences, headers = utils.load_fasta(ref_fasta)
    seq_lengths = {seq_id: len(seq) for seq_id, seq in sequences.items()}
    for seq_id, desc in headers:
        logger.info(f"  {seq_id}: {seq_lengths[seq_id]:,} bp")

    if circular:
        logger.info("  Circular genome mode: ON")

    # Stage 1: Detection
    if detection_input:
        logger.info("=" * 60)
        logger.info(f"Stage 1: Loading pre-computed coordinates from {detection_input}")
        logger.info("=" * 60)
        if not os.path.isfile(detection_input):
            logger.error(f"  Detection input file not found: {detection_input}")
            _save_results([], output_dir)
            return []
        tandems = _load_detection_input_tsv(
            detection_input, sequences, seq_lengths,
            circular=circular, flag_distance=flag_distance
        )
    else:
        # Decompress reference if gzipped (NUCmer can't read .gz)
        nucmer_ref, needs_cleanup = utils.decompress_if_gzipped(
            ref_fasta, output_dir=str(output_dir)
        )

        try:
            # Run NUCmer self-alignment
            coords_file = run_nucmer_self_alignment(
                nucmer_ref, str(output_dir), threads=threads
            )
        finally:
            # Clean up decompressed file
            if needs_cleanup and os.path.exists(nucmer_ref):
                os.unlink(nucmer_ref)

        # Parse alignments
        alignments = parse_coords(coords_file)
        logger.info(f"  Parsed {len(alignments)} alignments")

        # Identify tandem duplications (with sequences for complexity check)
        tandems = identify_tandem_duplications(
            alignments,
            seq_lengths=seq_lengths,
            max_distance=max_distance,
            min_identity=min_identity,
            min_size=min_size,
            circular=circular,
            sequences=sequences,
            apply_short_dup_filters=apply_short_dup_filters,
        )

    if not tandems:
        logger.info("  No tandem duplications detected.")
        _save_results([], output_dir)
        return []

    # Stage 2: Optional CDS subcategorization
    if annotation:
        logger.info("=" * 60)
        logger.info(f"Stage 2: CDS subcategorization from {annotation}")
        logger.info("=" * 60)
        if not os.path.isfile(annotation):
            logger.warning(
                f"  Annotation file not found: {annotation}. Skipping Stage 2."
            )
        else:
            tandems = subcategorize_by_annotation(tandems, annotation)
    else:
        logger.info(
            "Stage 2 skipped (no annotation provided; use --annotation to enable)"
        )

    # Stage 3: Classify mechanisms
    logger.info("=" * 60)
    logger.info("Stage 3: Mechanism classification")
    logger.info("=" * 60)
    classified = classify_tandem_mechanisms(
        tandems, sequences, flag_distance=flag_distance,
        hr_min_consec=hr_min_consec,
        hr_min_identity=hr_min_identity,
        hr_max_inward=hr_max_inward,
        hr_outward_ext=hr_outward_ext,
        hr_complexity_filter=hr_complexity_filter,
    )

    # Save results
    _save_results(classified, output_dir)

    return classified


def _save_results(tandems, output_dir):
    """Save module 1 results to TSV and JSON.

    Output columns are grouped by stage:
      Stage 1 (always): seq_id, copy coords, distance, size, identity
      Stage 2 (if annotation provided): n_cds_*, content_category, boundary_type
      Stage 3 (always): HR/microhomology classification + mechanism_confidence
    """
    tsv_path = os.path.join(output_dir, "tandem_duplications.tsv")
    with open(tsv_path, "w") as f:
        header = [
            # Stage 1 — Detection
            "seq_id", "copy1_start", "copy1_end", "copy2_start", "copy2_end",
            "distance", "nucmer_overlap_bp", "circular_closer",
            "size", "identity", "orientation", "proximity_class",
            "distance_to_size_ratio", "detection_source",
            # Stage 2 — CDS subcategorization (NA if no annotation)
            "n_cds_copy1", "n_cds_copy2", "n_cds_gap",
            "content_category", "boundary_type",
            # Stage 3 — Mechanism classification
            "is_hr_signature", "hr_match_len", "hr_identity", "hr_scenario",
            "microhomology_bp", "microhomology_seq",
            "junction_gap_bp", "junction_seq",
            "mechanism_confidence",
            "classification_note",
        ]
        f.write("\t".join(header) + "\n")
        for td in tandems:
            row = [str(td.get(h, "")) for h in header]
            f.write("\t".join(row) + "\n")

    logger.info(f"  Results saved to {tsv_path}")

    json_path = os.path.join(output_dir, "tandem_duplications.json")
    with open(json_path, "w") as f:
        json.dump(tandems, f, indent=2)

    logger.info(f"  JSON results saved to {json_path}")
