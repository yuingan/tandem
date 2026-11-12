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
        inner junction. Reports search_confidence (low/moderate/high)
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
from collections import Counter
from pathlib import Path

from . import utils
from .junction import KmerIndex, classify_junction, HR_MIN_CONSECUTIVE
from .repeats import annotate_repeat_families

logger = logging.getLogger("tandem")


# Mechanism confidence thresholds (copy length in bp)
CONFIDENCE_LOW_MAX = 300      # <300bp: flanking windows overlap copy itself
CONFIDENCE_MODERATE_MAX = 500 # 300-500bp: moderate
# >=500bp: high

# The two copies of a tandem duplication sit next to each other. A wider
# spacer means two dispersed copies of the same repeat, which is a substrate
# for a duplication rather than the product of one, and the sequence between
# them is not a junction. Some slack for filler bases and small indels.
MAX_TANDEM_GAP = 100

# Junction signatures are erased by substitution. Once the copies have drifted
# below this identity the microhomology window no longer reflects the event
# that made it, and the decay runs one way: HR signatures degrade into
# non-HR-looking ones, never the reverse.
DIVERGED_BELOW = 99.0

# A copy found at this many loci is a repeat family (IS element, rRNA operon)
# rather than a one-off duplication.
REPEAT_FAMILY_MIN_LOCI = 3


def _nucmer_has_threads():
    """MUMmer 4 accepts -t, MUMmer 3 does not."""
    import subprocess
    try:
        out = subprocess.run(
            ["nucmer", "--help"], capture_output=True, text=True, timeout=10
        )
    except (OSError, subprocess.SubprocessError):
        return False
    return "--threads" in (out.stdout + out.stderr)


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
        # MUMmer's docs require this for self-alignment. Without it MUMmer 3
        # drops every repeat as shadowed by the full-length diagonal. MUMmer 4
        # keeps most of them only because its shadow check never looks at
        # the first alignment, which is the diagonal; it still drops repeats
        # shadowed by any other alignment.
        "--nosimplify",
        "-l", str(min_match),
        "-c", str(min_cluster),
        "--prefix", prefix,
    ]
    # -t arrived in MUMmer 4. On MUMmer 3 it aborts with "Unknown option".
    if threads > 1 and _nucmer_has_threads():
        cmd[1:1] = ["-t", str(threads)]
    elif threads > 1:
        logger.warning(
            "  NUCmer does not accept -t; running single-threaded. "
            "Install MUMmer >= 4 for a parallel self-alignment."
        )
    cmd += [str(ref_fasta), str(ref_fasta)]
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


def _units(d, unit):
    """How many units make period d, give or take indels; 0 if no whole number does."""
    k = round(d / unit)
    if k >= 1 and abs(d - k * unit) <= min(max(10, 0.05 * d), unit / 4):
        return k
    return 0


def _array_regions(periodic, min_copies=2.5):
    """Tandem arrays from the self-alignments whose two sides overlap.

    Each entry (seq_id, start, end, period, identity) says start..end repeats
    with that period. A degenerate array seldom aligns to itself in one piece
    at its own period; the pieces at one, two, three... units cover it
    between them. Entries that overlap, at whole numbers of one unit, are
    merged into one region, shortest period first so the unit is the
    shortest that explains it. Pieces of one array split by up to two units
    of decayed sequence are joined. A region has to span min_copies units:
    below that it is two copies with a little more alignment past the
    junction, not an array.
    """
    regions = []
    for seq_id, start, end, p, identity in sorted(periodic, key=lambda e: e[3]):
        hits = [r for r in regions
                if r["seq_id"] == seq_id and _units(p, r["unit"])
                and start <= r["end"] + r["unit"] // 2
                and r["start"] - r["unit"] // 2 <= end]
        if not hits:
            regions.append({"seq_id": seq_id, "start": start, "end": end,
                            "unit": p, "identity": identity})
            continue
        keep = min(hits, key=lambda r: r["unit"])
        for r in hits:
            if r is not keep and _units(r["unit"], keep["unit"]):
                keep["start"] = min(keep["start"], r["start"])
                keep["end"] = max(keep["end"], r["end"])
                if _units(r["unit"], keep["unit"]) == 1:
                    keep["identity"] = max(keep["identity"], r["identity"])
                regions.remove(r)
        keep["start"] = min(keep["start"], start)
        keep["end"] = max(keep["end"], end)
        if _units(p, keep["unit"]) == 1:
            keep["identity"] = max(keep["identity"], identity)

    regions.sort(key=lambda r: (r["seq_id"], r["start"]))
    joined = []
    for r in regions:
        last = joined[-1] if joined else None
        if (last and last["seq_id"] == r["seq_id"]
                and _units(r["unit"], last["unit"]) == 1
                and r["start"] - last["end"] <= 2 * last["unit"]):
            last["end"] = max(last["end"], r["end"])
            last["identity"] = max(last["identity"], r["identity"])
        else:
            joined.append(r)
    for r in joined:
        r["copies"] = (r["end"] - r["start"] + 1) / r["unit"]
    return [r for r in joined if r["copies"] >= min_copies]


def _array_call(region):
    """One call for a whole array, its first two units standing in as copies."""
    start, unit = region["start"], region["unit"]
    copies = round(region["copies"], 1)
    if copies == int(copies):
        copies = int(copies)
    return {
        "seq_id": region["seq_id"],
        "copy1_start": start,
        "copy1_end": start + unit - 1,
        "copy2_start": start + unit,
        "copy2_end": start + 2 * unit - 1,
        "distance": 0,
        "nucmer_overlap_bp": 0,
        "circular_closer": False,
        "size": unit,
        "copy_number": copies,
        "identity": region["identity"],
        "orientation": "direct",
        "structure": "tandem",
        "distance_to_size_ratio": 0.0,
        "detection_source": "native",
    }


def _period(call):
    """Distance from the start of one copy to the start of the next."""
    return call["copy2_start"] - call["copy1_start"] - call["nucmer_overlap_bp"]


def _extent(call):
    """First and last base of a call, all copies included."""
    if call["copy_number"] > 2:
        return call["copy1_start"], call["copy1_start"] + call["size"] * call["copy_number"] - 1
    return call["copy1_start"], call["copy2_end"]


def _mark_unit_multiples(calls, regions):
    """Set aside pair calls that only restate an array.

    Inside an array of six 113 bp units, the self-alignment offset by 339 bp
    looks like a two-copy tandem, and a piece of the alignment at 113 bp
    like another, though neither is a duplication of its own. A pair lying
    inside an array region, at a whole number of its units, becomes a
    `unit_multiple` and records the unit. Regions whose unit is too short to
    report as an array still count.
    """
    by_seq = {}
    for r in sorted(regions, key=lambda r: r["unit"]):  # shortest unit first
        by_seq.setdefault(r["seq_id"], []).append(r)

    n = 0
    for call in calls:
        if call["structure"] not in ("tandem", "dispersed_pair"):
            continue
        d = _period(call)
        a, b = _extent(call)
        for r in by_seq.get(call["seq_id"], []):
            slack = r["unit"] // 2
            if r["start"] - slack <= a and b <= r["end"] + slack and _units(d, r["unit"]):
                call["structure"] = "unit_multiple"
                call["unit"] = r["unit"]
                n += 1
                break
    return n


def _sides(call):
    """The two stretches a call aligns to each other. For an array that is
    every unit but the last against every unit but the first."""
    if call["copy_number"] > 2:
        start, end = _extent(call)
        return (start, end - call["size"]), (start + call["size"], end)
    return ((call["copy1_start"], call["copy1_end"]),
            (call["copy2_start"], call["copy2_end"]))


def _covered(a, b):
    """Fraction of stretch a that lies inside stretch b."""
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1) / (a[1] - a[0] + 1)


def _drop_same_pair(calls, slack=20):
    """Keep one call per pair of copies.

    NUCmer aligns every pair of copies twice, once from each side. Once the
    copies have diverged the two alignments can stop at different bases, and
    the same pair turns up twice with slightly different ends. Two direct
    calls with the same period whose sides mostly overlap are one pair; the
    longer one stays. Copies of a family spaced evenly along the genome have
    the same period too, but their sides don't overlap.
    """
    direct = ("tandem", "dispersed_pair")
    kept = []
    for call in sorted(calls, key=lambda c: _extent(c)[0] - _extent(c)[1]):
        if call["structure"] in direct:
            d = _period(call)
            one, two = _sides(call)
            same = False
            for other in kept:
                if (other["seq_id"] != call["seq_id"]
                        or other["structure"] not in direct
                        or abs(_period(other) - d) > slack):
                    continue
                other_one, other_two = _sides(other)
                if _covered(one, other_one) >= 0.5 and _covered(two, other_two) >= 0.5:
                    same = True
                    break
            if same:
                continue
        kept.append(call)
    return kept, len(calls) - len(kept)


def identify_tandem_duplications(alignments, seq_lengths, max_distance=50000,
                                  min_identity=80.0, min_size=200,
                                  circular=True, sequences=None,
                                  short_dup_threshold=SHORT_DUP_THRESHOLD,
                                  apply_short_dup_filters=True,
                                  max_tandem_gap=MAX_TANDEM_GAP):
    """Identify tandem duplications from NUCmer self-alignment results.

    Every hit is labelled `tandem` or `dispersed_pair` depending on the gap
    between the copies, `inverted_pair` when one copy is reversed, or
    `unit_multiple` when its period is a whole multiple of a shorter repeat
    unit in the same region. Only tandems carry a junction; the rest are
    reported separately and are not duplications. A pair aligned twice is
    reported once.

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
    """
    tandems = []
    periodic = []
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

        inverted = (aln["ref_start"] < aln["ref_end"]) != (aln["qry_start"] < aln["qry_end"])

        # Calculate distance and assign copy order
        nucmer_overlap_bp = 0
        unit_size = aln["ref_alen"]

        if pos1[1] < pos2[0]:
            linear_distance = pos2[0] - pos1[1] - 1
            copy1, copy2 = pos1, pos2
        elif pos2[1] < pos1[0]:
            linear_distance = pos1[0] - pos2[1] - 1
            copy1, copy2 = pos2, pos1
        else:
            # The two sides overlap. A little is the microhomology at a
            # two-copy junction; more means the stretch repeats with this
            # period, and arrays are built from those after the loop.
            # Reversed and overlapping is one stretch folding back on itself.
            if inverted:
                continue
            overlap = min(pos1[1], pos2[1]) - max(pos1[0], pos2[0]) + 1
            first, second = (pos1, pos2) if pos1[0] <= pos2[0] else (pos2, pos1)
            if overlap > aln["ref_alen"] * 0.1:
                period = second[0] - first[0]
                # A repeat unit shifts both ends of the alignment by the same
                # amount. One that starts on the diagonal and drifts off it
                # through a low-complexity stretch says nothing about a unit.
                if period >= 10 and abs(second[1] - first[1] - period) <= max(10, period // 10):
                    periodic.append((aln["ref_id"], first[0], second[1], period, aln["identity"]))
                continue
            nucmer_overlap_bp = overlap
            copy1, copy2 = list(first), list(second)
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

        orientation = "inverted" if inverted else "direct"

        # decides whether the sequence between the copies is a junction at
        # all. Reversed copies have none however close they are: the end of
        # one never meets the start of the other.
        if inverted:
            structure = "inverted_pair"
        elif nucmer_overlap_bp > 0 or distance <= max_tandem_gap:
            structure = "tandem"
        else:
            structure = "dispersed_pair"

        # Deduplicate
        key = (seq_id, copy1[0], copy1[1], copy2[0], copy2[1])
        if key in seen:
            continue
        seen.add(key)

        # Short-duplication filters (<short_dup_threshold bp)
        dup_size = unit_size
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
            "size": unit_size,
            "copy_number": 2,
            "identity": aln["identity"],
            "orientation": orientation,
            "structure": structure,
            "distance_to_size_ratio": round(distance / unit_size, 2) if unit_size > 0 else 0.0,
            "detection_source": "native",
        })

    regions = _array_regions(periodic)
    tandems, n_twice = _drop_same_pair(tandems)
    n_multiples = _mark_unit_multiples(tandems, regions)

    # Arrays whose unit is below the size floor aren't duplications to report,
    # but they are listed so they can be counted.
    n_arrays = n_short = 0
    for r in regions:
        if apply_short_dup_filters and r["unit"] < short_dup_threshold and sequences is not None:
            ref = sequences.get(r["seq_id"], "")
            if ref and not _sequence_complexity_ok(ref[r["start"] - 1:r["start"] - 1 + r["unit"]]):
                n_filtered_by_complexity += 1
                continue
        call = _array_call(r)
        if r["unit"] < min_size:
            call["structure"] = "short_unit_array"
            n_short += 1
        else:
            n_arrays += 1
        tandems.append(call)

    if n_twice:
        logger.info(f"  {n_twice} call(s) dropped as the same pair aligned twice")
    if n_arrays or n_short:
        logger.info(f"  {n_arrays} tandem array(s), and {n_short} with a unit below {min_size}bp")
    if n_multiples:
        logger.info(f"  {n_multiples} call(s) are whole multiples of an array's unit")

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

    n = Counter(t["structure"] for t in tandems)
    logger.info(
        f"  {n['tandem']} tandem duplication candidate(s), "
        f"{n['dispersed_pair']} dispersed repeat pair(s) (gap > {max_tandem_gap}bp), "
        f"{n['inverted_pair']} inverted pair(s), "
        f"{n['unit_multiple']} multiple(s) of a shorter unit, "
        f"{n['short_unit_array']} array(s) with a unit below {min_size}bp"
    )
    return tandems


def describe_junctions(tandems, sequences, flag_distance=2000,
                       hr_min_consec=None, hr_min_identity=None,
                       hr_max_inward=None, hr_outward_ext=None,
                       hr_complexity_filter=True,
                       diverged_below=DIVERGED_BELOW):
    """Describe the sequence architecture at each junction.

    Records whether an extended direct repeat flanks the duplicated segment
    (and in which geometry), and the microhomology at the inner junction.

    These are observations about substrate, not mechanism calls. What a
    reference genome shows is what is there now: substitution rewrites
    junctions, and gene conversion between tandem copies can hold a flanking
    repeat identical long after the event. Copies below diverged_below are
    flagged, not filtered.

    Notes column: junction_distant, too_small_for_repeat_search, copies_diverged
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
    tally = {"flanking_repeat": 0, "no_flanking_repeat": 0,
             "array": 0, "dispersed_pair": 0}
    n_diverged = 0

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

        seen = classify_junction(seq, c1s_0, c1e_0, c2s_0, c2e_0,
                                 overlap=td.get("nucmer_overlap_bp", 0), **hr_kwargs)
        for key in ("flanking_repeat", "repeat_len", "repeat_identity",
                    "repeat_geometry", "microhomology_bp", "microhomology_seq",
                    "junction_gap_bp", "junction_seq"):
            td[key] = seen[key]

        # how far out the flanking search could look; not a biological claim
        copy_len = td["size"]
        if copy_len < CONFIDENCE_LOW_MAX:
            td["search_confidence"] = "low"
        elif copy_len < CONFIDENCE_MODERATE_MAX:
            td["search_confidence"] = "moderate"
        else:
            td["search_confidence"] = "high"

        # Arrays and dispersed pairs get their own labels: an array is
        # flanked by a full copy of itself at every internal junction, and a
        # dispersed pair has no junction.
        if td.get("structure", "tandem") != "tandem":
            td["architecture"] = "dispersed_pair"
        elif td.get("copy_number", 2) > 2:
            td["architecture"] = "array"
        elif td["flanking_repeat"]:
            td["architecture"] = "flanking_repeat"
        else:
            td["architecture"] = "no_flanking_repeat"

        td["copies_diverged"] = td["identity"] < diverged_below

        notes = []
        if td["copies_diverged"]:
            notes.append("copies_diverged")
        if td["distance"] > flag_distance:
            notes.append("junction_distant")
        if td["size"] < 2 * HR_MIN_CONSECUTIVE:
            notes.append("too_small_for_repeat_search")
        td["classification_note"] = ";".join(notes)

        classified.append(td)
        tally[td["architecture"]] += 1
        if td["copies_diverged"] and td["architecture"] in ("flanking_repeat", "no_flanking_repeat"):
            n_diverged += 1

    logger.info(
        f"  Junction architecture: {tally['flanking_repeat']} with flanking "
        f"direct repeat, {tally['no_flanking_repeat']} without, "
        f"{tally['array']} array(s), {tally['dispersed_pair']} dispersed pair(s)"
    )
    if n_diverged:
        logger.info(
            f"  {n_diverged} of these have copies below {diverged_below}% identity; "
            f"the repeat observation describes the present substrate, not formation"
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
        if not td["flanking_repeat"]:
            mh = td["microhomology_bp"]
            mh_dist[mh] = mh_dist.get(mh, 0) + 1
    if mh_dist:
        top_mh = sorted(mh_dist.items(), key=lambda x: -x[1])[:10]
        logger.info(f"  Non-HR microhomology distribution (top 10):")
        for mh_bp, count in top_mh:
            logger.info(f"    {mh_bp} bp: {count}")

    return classified


def _parse_gff_cds_features(annotation_path):
    """Parse a GFF3 file and return a list of CDS features."""
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
    """Parse a GenBank file and return a list of CDS features."""
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
    """Count CDS fully contained (or >= min_overlap_frac overlapping) in interval."""
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
    """Determine whether a boundary falls inside a CDS, at a CDS edge, or in IGS."""
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
                hr_complexity_filter=True,
                max_tandem_gap=MAX_TANDEM_GAP,
                diverged_below=DIVERGED_BELOW,
                repeat_family_min_loci=REPEAT_FAMILY_MIN_LOCI,
                keep_dispersed=False):
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
            max_tandem_gap=max_tandem_gap,
        )
        annotate_repeat_families(
            tandems, alignments, min_loci=repeat_family_min_loci
        )

    # Not duplications; mixing them into the main table is what makes an IS
    # pair a few kb apart look like a result. Both files are written even when
    # empty, so a rerun into the same folder can't leave an old one behind.
    dispersed = [t for t in tandems if t.get("structure") == "dispersed_pair"]
    inverted = [t for t in tandems if t.get("structure") == "inverted_pair"]
    multiples = [t for t in tandems if t.get("structure") == "unit_multiple"]
    short_arrays = [t for t in tandems if t.get("structure") == "short_unit_array"]
    _save_pairs(dispersed, output_dir, "dispersed_pairs.tsv")
    _save_pairs(inverted, output_dir, "inverted_pairs.tsv")
    _save_pairs(multiples, output_dir, "unit_multiples.tsv", extra=("unit",))
    _save_pairs(short_arrays, output_dir, "short_unit_arrays.tsv", extra=("copy_number",))
    if dispersed and not keep_dispersed:
        logger.info(
            f"  {len(dispersed)} dispersed pair(s) written to "
            f"dispersed_pairs.tsv (use --keep-dispersed to merge them in)"
        )
    if inverted:
        logger.info(f"  {len(inverted)} inverted pair(s) written to inverted_pairs.tsv")
    if multiples:
        logger.info(f"  {len(multiples)} multiple(s) of a shorter unit written to unit_multiples.tsv")
    if short_arrays:
        logger.info(f"  {len(short_arrays)} array(s) with a unit below {min_size}bp written "
                    f"to short_unit_arrays.tsv")

    keep = {"tandem", None}
    if keep_dispersed:
        keep.add("dispersed_pair")
    tandems = [t for t in tandems if t.get("structure") in keep]

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
    classified = describe_junctions(
        tandems, sequences, flag_distance=flag_distance,
        hr_min_consec=hr_min_consec,
        hr_min_identity=hr_min_identity,
        hr_max_inward=hr_max_inward,
        hr_outward_ext=hr_outward_ext,
        hr_complexity_filter=hr_complexity_filter,
        diverged_below=diverged_below,
    )

    # Save results
    _save_results(classified, output_dir)

    return classified


def _save_pairs(pairs, output_dir, name, extra=()):
    """Repeat pairs that are not tandem duplications.

    Direct copies too far apart are the substrate an unequal exchange acts
    on; inverted copies are what an inversion needs. Neither is a duplication,
    and nor is a multiple of a shorter unit inside an array. Arrays whose unit
    is below the size floor go through here too.
    """
    path = os.path.join(output_dir, name)
    header = [
        "seq_id", "copy1_start", "copy1_end", "copy2_start", "copy2_end",
        "distance", "size", "identity", "orientation",
        "repeat_loci", "dispersed_family", *extra,
    ]
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for pair in pairs:
            f.write("\t".join(str(pair.get(h, "")) for h in header) + "\n")


def _save_results(tandems, output_dir):
    """Save module 1 results to TSV and JSON.

    Output columns are grouped by stage:
      Stage 1 (always): seq_id, copy coords, distance, size, identity
      Stage 2 (if annotation provided): n_cds_*, content_category, boundary_type
      Stage 3 (always): HR/microhomology classification + search_confidence
    """
    tsv_path = os.path.join(output_dir, "tandem_duplications.tsv")
    with open(tsv_path, "w") as f:
        header = [
            # Stage 1 — Detection
            "seq_id", "copy1_start", "copy1_end", "copy2_start", "copy2_end",
            "distance", "nucmer_overlap_bp", "circular_closer",
            "size", "copy_number", "identity", "orientation", "structure",
            "distance_to_size_ratio", "repeat_loci", "dispersed_family",
            "detection_source",
            # Stage 2 — CDS subcategorization (NA if no annotation)
            "n_cds_copy1", "n_cds_copy2", "n_cds_gap",
            "content_category", "boundary_type",
            # Stage 3 — Mechanism classification
            "flanking_repeat", "repeat_len", "repeat_identity", "repeat_geometry",
            "microhomology_bp", "microhomology_seq",
            "junction_gap_bp", "junction_seq",
            "architecture", "copies_diverged", "search_confidence",
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
