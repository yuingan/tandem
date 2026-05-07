"""Core junction analysis logic for Tandem.

This module provides:
- HR detection via flanking repeat alignment (matching Part 3 approach)
- Microhomology measurement at junctions (reports bp count, no binning)
- Junction candidate generation for modules 2 & 3
- Extend-to-unique algorithm for creating diagnostic junction references
- K-mer indexing for fast uniqueness checking
"""

import logging

import numpy as np

logger = logging.getLogger("tandem")


# =============================================================================
# HR detection parameters (matching Part 3)
# =============================================================================

HR_MIN_CONSECUTIVE = 35     # bp - minimum consecutive match length
HR_MIN_IDENTITY = 0.92      # minimum alignment identity
HR_MAX_WINDOW = 2000        # bp - cap for adaptive window (inside extent)
HR_OUTWARD_EXTENSION = 200  # bp - extend windows OUTWARD beyond copy boundaries
                            # to tolerate NUCmer-reported boundary uncertainty.
                            # The R region mediating HR may start a few bp
                            # inside/outside the NUCmer boundary.

# Junction extraction
JUNCTION_FLANK = 500        # bp flanking sequence for microhomology check

# HR match complexity filter
HR_COMPLEXITY_K = 4         # k-mer size for complexity check
HR_COMPLEXITY_MIN = 0.20    # minimum fraction of possible k-mers


def _hr_sequence_complexity_ok(seq, k=HR_COMPLEXITY_K,
                                min_fraction=HR_COMPLEXITY_MIN):
    """Check if sequence has enough k-mer diversity to be an HR mediator.

    Simple repeats like (GGCCTTAA)ₙ are tandem microsatellites that
    can mediate replication slippage, not HR. A real HR-mediating R
    should be complex enough to be a unique genomic locus.

    Returns True if sequence has adequate complexity.
    """
    if len(seq) < k:
        return True  # too short to assess
    kmers = set()
    for i in range(len(seq) - k + 1):
        kmers.add(seq[i:i + k])
    max_possible = min(4 ** k, len(seq) - k + 1)
    return (len(kmers) / max_possible) >= min_fraction if max_possible > 0 else True

# Extend-to-unique parameters (modules 2 & 3)
EXTEND_START_K = 15         # initial extension on each side
EXTEND_MAX_K = 75           # default max extension (2*75=150bp fits Illumina 150bp reads)


# =============================================================================
# K-mer index for fast uniqueness checking
# =============================================================================

class KmerIndex:
    """Index of k-mers in a reference genome for fast uniqueness checking."""

    def __init__(self, sequences, min_k=EXTEND_START_K):
        self.sequences = sequences
        self.min_k = min_k
        self._concat_seq = ""
        self._seq_offsets = {}
        offset = 0
        for seq_id, seq in sequences.items():
            self._seq_offsets[seq_id] = offset
            self._concat_seq += seq + "N" * 100
            offset += len(seq) + 100

        logger.debug(
            f"KmerIndex: {len(sequences)} sequences, "
            f"{len(self._concat_seq):,} bp total"
        )

    def count_occurrences(self, query, max_count=None):
        """Count exact occurrences in the reference (both strands).

        Args:
            query: DNA sequence string
            max_count: stop counting after this many (for early exit)

        Returns:
            Number of exact occurrences (or max_count if reached early)
        """
        count = 0
        start = 0
        seq = self._concat_seq
        q = query.upper()

        while True:
            pos = seq.find(q, start)
            if pos == -1:
                break
            count += 1
            if max_count is not None and count >= max_count:
                return count
            start = pos + 1

        rc = reverse_complement(q)
        if rc != q:
            start = 0
            while True:
                pos = seq.find(rc, start)
                if pos == -1:
                    break
                count += 1
                if max_count is not None and count >= max_count:
                    return count
                start = pos + 1

        return count

    def is_unique(self, query):
        """Check if query occurs at most once (stops after finding 2)."""
        return self.count_occurrences(query, max_count=2) <= 1

    def is_absent(self, query):
        """Check if query has zero occurrences (stops after finding 1)."""
        return self.count_occurrences(query, max_count=1) == 0


def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


# =============================================================================
# HR detection (matching Part 3 approach)
# =============================================================================

def _get_adaptive_window(dup_length):
    """Adaptive window sizing based on duplication length.

    Cap at dup_length // 2 so at least half the comparison is
    genuinely flanking sequence, not the copy itself.
    Large dups (>=2kb): window = 2kb cap.
    """
    if dup_length <= 0:
        return HR_MAX_WINDOW
    # Use at most half the copy length, so the check region
    # extends outside the copy
    half_dup = dup_length // 2
    return min(half_dup, HR_MAX_WINDOW)


def _parse_alignment_format(alignment):
    """Parse Bio.Align.PairwiseAligner format string into aligned sequences.

    Biopython's format() wraps every 60 characters. The format varies
    slightly by version:

    Most blocks (3 tokens):
        target  0 ACGTACGT...
                0 ||||||||...
        query   0 ACGTACGT...

    Last block (4 tokens, includes end position):
        target 180 ACGTACGT... 192
               180 ||||||||... 192
        query  180 ACGTACGT... 192

    In both cases, the sequence is at index 2 (parts[2]).

    Returns: (target_aligned, query_aligned) as strings with gaps
    """
    lines = alignment.format().strip().split("\n")

    target_parts = []
    query_parts = []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue

        parts = stripped.split()
        if len(parts) < 3:
            continue

        # Target lines start with 'target', query lines start with 'query'
        # The sequence is always at index 2
        if parts[0] == "target":
            target_parts.append(parts[2])
        elif parts[0] == "query":
            query_parts.append(parts[2])
        # Match lines (||||) and other lines are skipped

    return "".join(target_parts), "".join(query_parts)


def _calculate_alignment_metrics(target_aln, query_aln):
    """Calculate metrics from aligned sequence strings.
    Returns: (total_matches, identity, max_consecutive_matches)
    """
    total_matches = 0
    aligned_length = 0
    max_consecutive = 0
    current_consecutive = 0

    for a, b in zip(target_aln, query_aln):
        if a != '-' or b != '-':
            aligned_length += 1

        if a == b and a != '-':
            total_matches += 1
            current_consecutive += 1
            max_consecutive = max(max_consecutive, current_consecutive)
        else:
            current_consecutive = 0

    identity = total_matches / aligned_length if aligned_length > 0 else 0.0
    return total_matches, identity, max_consecutive


def _best_identity_window(target_aln, query_aln, min_window=35,
                           min_identity=0.92):
    """Find the longest sub-region of the alignment with identity >= threshold.

    Slides a window of `min_window` across the alignment, checking local
    identity. If a passing window is found, it is extended in both
    directions as long as the extended region still meets the threshold.

    This is the Option 2 (windowed identity / HSP-style) approach:
    it correctly detects diverged HR repeats where mismatches are scattered
    (breaking consecutive-match runs) but local identity remains high.

    Args:
        target_aln: aligned target string (may contain '-' gaps)
        query_aln: aligned query string (may contain '-' gaps)
        min_window: minimum window size to test (default 35)
        min_identity: minimum identity threshold within the window

    Returns:
        (best_length, best_identity) — length and identity of the best
        passing window. Returns (0, 0.0) if no window meets the threshold.
    """
    n = len(target_aln)
    if n < min_window:
        return 0, 0.0

    # Build per-column match array (1 = match, 0 = mismatch or gap)
    # and aligned array (1 = aligned column, 0 = double-gap which shouldn't happen)
    match_arr = []
    for a, b in zip(target_aln, query_aln):
        match_arr.append(1 if (a == b and a != '-') else 0)

    # Prefix sums for O(1) window queries
    prefix = [0] * (n + 1)
    for i in range(n):
        prefix[i + 1] = prefix[i] + match_arr[i]

    def _window_identity(start, end):
        """Identity of alignment columns [start, end)."""
        length = end - start
        if length <= 0:
            return 0, 0.0
        matches = prefix[end] - prefix[start]
        return length, matches / length

    # Phase 1: Find any passing window of exactly min_window size
    best_start = -1
    best_ident = 0.0
    for start in range(n - min_window + 1):
        end = start + min_window
        _, ident = _window_identity(start, end)
        if ident >= min_identity and ident > best_ident:
            best_ident = ident
            best_start = start

    if best_start < 0:
        return 0, 0.0

    # Phase 2: Extend the best window in both directions while identity holds
    ext_start = best_start
    ext_end = best_start + min_window

    # Extend rightward
    while ext_end < n:
        _, new_ident = _window_identity(ext_start, ext_end + 1)
        if new_ident >= min_identity:
            ext_end += 1
        else:
            break

    # Extend leftward
    while ext_start > 0:
        _, new_ident = _window_identity(ext_start - 1, ext_end)
        if new_ident >= min_identity:
            ext_start -= 1
        else:
            break

    final_length = ext_end - ext_start
    _, final_identity = _window_identity(ext_start, ext_end)
    return final_length, final_identity


def _run_local_alignment(seq1, seq2, min_window=HR_MIN_CONSECUTIVE,
                         min_identity=HR_MIN_IDENTITY):
    """Run local alignment and compute windowed identity metrics.

    Uses strict scoring (match=+2, mismatch=-3, gap_open=-5, gap_ext=-2)
    to constrain the alignment to the homology core, then applies
    sliding-window identity analysis (Option 2 / HSP-style) to find
    the best sub-region meeting the identity threshold.

    Returns: (window_length, window_identity, max_consecutive_matches)
        window_length: length of the best passing identity window (0 if none)
        window_identity: identity within that window (0.0 if none)
        max_consecutive_matches: longest uninterrupted exact match run
                                 (kept for backward compatibility / logging)
    """
    if len(seq1) == 0 or len(seq2) == 0:
        return 0, 0.0, 0

    try:
        from Bio.Align import PairwiseAligner
        aligner = PairwiseAligner()
        aligner.mode = 'local'
        aligner.match_score = 2
        aligner.mismatch_score = -3
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -2
        alignments = aligner.align(seq1, seq2)
        if alignments:
            target_aln, query_aln = _parse_alignment_format(alignments[0])
            if target_aln and query_aln:
                # Windowed identity (primary metric for HR detection)
                win_len, win_ident = _best_identity_window(
                    target_aln, query_aln,
                    min_window=min_window,
                    min_identity=min_identity,
                )
                # Max consecutive (secondary metric, kept for logging)
                _, _, max_consec = _calculate_alignment_metrics(
                    target_aln, query_aln
                )
                return win_len, win_ident, max_consec
        return 0, 0.0, 0
    except Exception as e:
        logger.debug(f"Alignment failed: {e}")
        return 0, 0.0, 0


def check_hr(genome_seq, c1s_0, c1e_0, c2s_0, c2e_0,
             min_consec=HR_MIN_CONSECUTIVE, min_identity=HR_MIN_IDENTITY,
             max_inward=HR_MAX_WINDOW, outward_ext=HR_OUTWARD_EXTENSION,
             complexity_filter=True):
    """Detect HR via flanking repeat alignment (adaptive window).

    In an R-B-R-B-R structure, NUCmer reports two copies with ambiguous
    bracketing. We check four combinations:

    Scenario 1 — R at start of copies, HR repeat AFTER copy 2:
      [R₁-B₁]-[R₂-B₂]-R₃
      Check 1a: start of copy 1 (R₁) vs after copy 2 (R₃)
      Check 1b: start of copy 2 (R₂) vs after copy 2 (R₃)

    Scenario 2 — R at end of copies, HR repeat BEFORE copy 1:
      R₁-[B₁-R₂]-[B₂-R₃]
      Check 2a: end of copy 1 (R₂) vs before copy 1 (R₁)
      Check 2b: end of copy 2 (R₃) vs before copy 1 (R₁)

    Window geometry — each window has a SEARCH direction (where R is
    expected, extends by `inward` bp) and a TOLERANCE direction (for
    NUCmer boundary uncertainty):

      OUTER boundaries where extension goes OUTSIDE the duplication
      (c1s_0 leftward, c2e_0 rightward): tolerance = outward_ext.
      User can increase this freely since it doesn't enter any copy body.

      OUTER boundaries where extension goes INTO a copy body
      (win_after_c2 leftward into copy 2, win_before_c1 rightward into
      copy 1): tolerance capped at MAX_BODY_TOLERANCE (200bp) regardless
      of --hr-outward-ext. Without this cap, large outward_ext would
      extend deep into body, and since body₁ ≡ body₂ in a tandem, the
      aligner would match body-tail vs body-tail → false positive HR.

      INNER boundaries (c2s_0 for check 1b, c1e_0 for check 2a): extend
      ONLY inward into the copy, zero extension across the inner junction.

    Threshold (windowed identity / HSP-style):
      Within each pairwise alignment, find the longest sub-region of
      length >= min_consec where identity >= min_identity.

    Guard: Skip if copy is too small (< 2 * min_consec).

    Returns: (is_hr, match_length, identity, scenario)
    """
    MAX_BODY_TOLERANCE = 200  # Never extend deeper into a copy body

    dup_length = c1e_0 - c1s_0

    if dup_length < 2 * min_consec:
        return False, 0, 0.0, "too_small_for_HR_check"

    inward = min(dup_length // 2, max_inward)
    seq_len = len(genome_seq)

    # Cap into-body tolerance: prevents body₁≡body₂ false matches
    body_tol = min(outward_ext, MAX_BODY_TOLERANCE)

    def _window(left, right):
        """Extract genome subsequence [left, right), clamped to genome bounds."""
        return genome_seq[max(0, left):min(seq_len, right)]

    def _check_pair(seq_a, seq_b, label):
        """Run alignment with windowed identity check + complexity filter."""
        if len(seq_a) < min_consec or len(seq_b) < min_consec:
            return None
        win_len, win_ident, _ = _run_local_alignment(
            seq_a, seq_b,
            min_window=min_consec,
            min_identity=min_identity,
        )
        if win_len >= min_consec and win_ident >= min_identity:
            # Complexity check: reject if BOTH windows are low-complexity.
            # A periodic body (e.g. GGCCTTAA×40) produces real alignments
            # from any two windows, but these are microsatellite periodicity,
            # not HR-mediating repeats. At least one window must contain
            # complex sequence for the match to be biologically meaningful.
            # Disabled by --no-hr-complexity-filter.
            if complexity_filter:
                if (not _hr_sequence_complexity_ok(seq_a) and
                        not _hr_sequence_complexity_ok(seq_b)):
                    return None
            return (True, win_len, win_ident,
                    f"{label}_w{inward}_ext{outward_ext}")
        return None

    # === Scenario 1: R at start of copies, R₃ is AFTER copy 2 ===

    # Window at c1s_0 (OUTER boundary, extension goes OUTSIDE):
    #   LEFT = outward_ext (outside dup, user-adjustable, no body issue)
    #   RIGHT = inward (into copy 1, searching for R₁)
    win_c1_start = _window(c1s_0 - outward_ext, c1s_0 + inward)

    # Window at c2s_0 (INNER boundary): ONLY extend right into copy 2
    win_c2_start = _window(c2s_0, c2s_0 + inward)

    # Window after c2e_0 (OUTER boundary, extension into body₂):
    #   LEFT = body_tol (into copy 2 body, CAPPED to prevent false match)
    #   RIGHT = inward (past copy 2, searching for R₃)
    win_after_c2 = _window(c2e_0 - body_tol, c2e_0 + inward)

    # Check 1a: R₁ (at outer 5' of copy 1) vs R₃ (after copy 2)
    result = _check_pair(win_c1_start, win_after_c2, "s1a_R1_vs_R3")
    if result:
        return result

    # Check 1b: R₂ (at inner start of copy 2) vs R₃ (after copy 2)
    result = _check_pair(win_c2_start, win_after_c2, "s1b_R2_vs_R3")
    if result:
        return result

    # === Scenario 2: R at end of copies, R₁ is BEFORE copy 1 ===

    # Window before c1s_0 (OUTER boundary, extension into body₁):
    #   LEFT = inward (before copy 1, searching for R₁)
    #   RIGHT = body_tol (into copy 1 body, CAPPED to prevent false match)
    win_before_c1 = _window(c1s_0 - inward, c1s_0 + body_tol)

    # Window at c1e_0 (INNER boundary): ONLY extend left into copy 1
    win_c1_end = _window(c1e_0 - inward, c1e_0)

    # Window at c2e_0 (OUTER boundary, extension goes OUTSIDE):
    #   LEFT = inward (into copy 2, searching for R₃)
    #   RIGHT = outward_ext (outside dup, user-adjustable, no body issue)
    win_c2_end = _window(c2e_0 - inward, c2e_0 + outward_ext)

    # Check 2a: R₂ (at inner end of copy 1) vs R₁ (before copy 1)
    result = _check_pair(win_c1_end, win_before_c1, "s2a_R2_vs_R1")
    if result:
        return result

    # Check 2b: R₃ (at outer 3' of copy 2) vs R₁ (before copy 1)
    result = _check_pair(win_c2_end, win_before_c1, "s2b_R3_vs_R1")
    if result:
        return result

    return False, 0, 0.0, "no_HR"


# =============================================================================
# Microhomology detection at junction
# =============================================================================

def find_microhomology(seq1, seq2, max_len=None):
    """Find maximum exact microhomology at junction.
    Checks if end of seq1 matches start of seq2.

    Args:
        seq1: upstream flank (end of copy1 region)
        seq2: downstream flank (start of copy2 region)
        max_len: maximum microhomology to check (default: full input length).
                 Previously defaulted to 60, which silently returned 0 for
                 true MH > 60bp because partial k-mer comparisons misalign.

    Returns:
        int: number of bp of microhomology
    """
    if max_len is None:
        max_check = min(len(seq1), len(seq2))
    else:
        max_check = min(max_len, len(seq1), len(seq2))
    for k in range(max_check, 0, -1):
        if seq1[-k:] == seq2[:k]:
            return k
    return 0


def classify_junction(genome_seq, c1s_0, c1e_0, c2s_0, c2e_0,
                      hr_min_consec=HR_MIN_CONSECUTIVE,
                      hr_min_identity=HR_MIN_IDENTITY,
                      hr_max_inward=HR_MAX_WINDOW,
                      hr_outward_ext=HR_OUTWARD_EXTENSION,
                      hr_complexity_filter=True):
    """Classify a tandem duplication junction.

    Two-step approach:
    1. Check for HR (flanking repeat alignment)
    2. Measure microhomology bp at junction

    Reports HR status and microhomology bp count.
    Does NOT bin into NHEJ/MMEJ/SSA categories.

    Args:
        genome_seq: contig sequence string
        c1s_0: 0-based start of copy1
        c1e_0: 0-based exclusive end of copy1
        c2s_0: 0-based start of copy2
        c2e_0: 0-based exclusive end of copy2
        hr_min_consec: minimum consecutive match length for HR detection
        hr_min_identity: minimum alignment identity for HR detection
        hr_max_inward: cap for inward window extent
        hr_outward_ext: bp to extend outward beyond copy boundary
        hr_complexity_filter: if True, reject HR calls where both windows
            are low-complexity (microsatellite periodicity). Default True.

    Returns:
        dict with classification details
    """
    # Step 1: Check HR
    is_hr, hr_match_len, hr_identity, hr_scenario = check_hr(
        genome_seq, c1s_0, c1e_0, c2s_0, c2e_0,
        min_consec=hr_min_consec,
        min_identity=hr_min_identity,
        max_inward=hr_max_inward,
        outward_ext=hr_outward_ext,
        complexity_filter=hr_complexity_filter,
    )

    # Step 2: Measure microhomology at junction
    flank_size = min(JUNCTION_FLANK, c1e_0 - c1s_0)
    upstream_flank = genome_seq[max(0, c1e_0 - flank_size) : c1e_0]
    downstream_flank = genome_seq[c2s_0 : min(len(genome_seq), c2s_0 + flank_size)]

    microhomology_bp = find_microhomology(upstream_flank, downstream_flank)
    microhomology_seq = upstream_flank[-microhomology_bp:] if microhomology_bp > 0 else ""

    # Junction gap (inter-copy distance)
    junction_gap = c2s_0 - c1e_0
    junction_seq = genome_seq[c1e_0:c2s_0] if junction_gap > 0 else ""

    return {
        "is_hr_signature": is_hr,
        "hr_match_len": hr_match_len,
        "hr_identity": round(hr_identity, 4),
        "hr_scenario": hr_scenario,
        "microhomology_bp": microhomology_bp,
        "microhomology_seq": microhomology_seq,
        "junction_gap_bp": max(0, junction_gap),
        "junction_seq": junction_seq[:100] if junction_seq else "",
    }


def classify_single_copy_junction(ref_seq, start_0, end_0,
                                   hr_min_consec=HR_MIN_CONSECUTIVE,
                                   hr_min_identity=HR_MIN_IDENTITY,
                                   hr_max_inward=HR_MAX_WINDOW,
                                   hr_outward_ext=HR_OUTWARD_EXTENSION,
                                   hr_complexity_filter=True):
    """Classify a tandem duplication junction from a SINGLE-COPY reference.

    For Modules 2 and 3: the reference genome contains only one copy of
    the duplicated region [start_0, end_0). The second copy exists only
    in the mutant/isolate.

    This function differs from classify_junction (used by Module 1) in
    two critical ways:

    Microhomology: The biological junction in the mutant is
      ref[end-k:end] + ref[start:start+k]  (end of Y followed by start
      of Y). So microhomology = overlap between END of Y and START of Y.
      classify_junction would incorrectly compare end of Y vs downstream
      genomic context.

    HR detection: Asks "does the reference have R-Y-R structure?" — are
      there direct repeats flanking Y on both sides? Compares the region
      around Y_start vs the region around Y_end. If they share homology,
      direct repeats flank Y, consistent with HR-mediated duplication.
      classify_junction would project a virtual copy 2 into non-homologous
      downstream sequence.

    Args:
        ref_seq: reference sequence string
        start_0: 0-based start of duplicated region
        end_0: 0-based exclusive end of duplicated region
        hr_min_consec: minimum window size for HR identity check
        hr_min_identity: minimum identity for HR detection
        hr_max_inward: cap for inward window extent
        hr_outward_ext: outward extension for NUCmer boundary tolerance

    Returns:
        dict with same keys as classify_junction for API compatibility
    """
    dup_length = end_0 - start_0
    seq_len = len(ref_seq)

    # --- Step 1: Microhomology at the biological junction ---
    # Junction in the mutant: ...ref[end-k:end] | ref[start:start+k]...
    # MH = how many bases at end of Y match start of Y
    flank_size = min(JUNCTION_FLANK, dup_length)
    upstream_flank = ref_seq[max(0, end_0 - flank_size) : end_0]
    downstream_flank = ref_seq[start_0 : min(seq_len, start_0 + flank_size)]

    microhomology_bp = find_microhomology(upstream_flank, downstream_flank)
    microhomology_seq = upstream_flank[-microhomology_bp:] if microhomology_bp > 0 else ""

    # --- Step 2: HR detection ---
    # Four scenarios test for direct repeats that could mediate HR.
    # All four search for a homology block ≥ hr_min_consec bp at
    # ≥ hr_min_identity identity.
    #
    # Scenario 3 (R-Y-R, primary): repeat OUTSIDE Y on both sides.
    #   Windows: [start-out, start+in] vs [end-in, end+out]
    #   Classical HR: rRNA operons or IS elements flanking a duplicated region.
    #
    # Scenario 1 (asymmetric): repeat OUTSIDE Y at start, INSIDE Y at end.
    #   Windows: [start-out, start+small_in] vs [end-in, end+small_in]
    #   Geometry: the duplication boundary at the end is INTERNAL to a
    #   repeat, with the second repeat copy upstream of Y_start in reference.
    #
    # Scenario 2 (asymmetric): repeat INSIDE Y at start, OUTSIDE Y at end.
    #   Windows: [start-small_in, start+in] vs [end-small_in, end+out]
    #   Geometry: the duplication boundary at the start is INTERNAL to a
    #   repeat, with the second repeat copy downstream of Y_end in reference.
    #
    # Scenario 4 (internal-internal): repeat INSIDE Y at both boundaries.
    #   Windows: [start, start+in] vs [end-in, end]
    #   Geometry: two homologous tracts (e.g. rRNA operons) sit at the
    #   internal start and end of Y. HR between them via sister-chromatid
    #   unequal crossing-over produces a tandem duplication of Y. This is
    #   distinct from the immediate-junction microhomology measurement —
    #   s4 looks for an EXTENDED repeat (≥ 35 bp at ≥ 92% identity) using
    #   the same windowed-HSP analysis as the other scenarios. The two
    #   windows must be disjoint, so inward is clamped to dup_length // 2.
    #
    # Adaptive outward extension: for large duplications (> 50kb), the
    # flanking R repeat may sit several kb outside the boundary (e.g.
    # rRNA operons spaced 5-10kb apart). The user-supplied hr_outward_ext
    # (default 200bp) is sufficient for boundary uncertainty but too small
    # for finding distant flanking repeats. Adaptive scaling: extend outward
    # up to 5kb for duplications ≥ 50kb.
    is_hr = False
    hr_match_len = 0
    hr_identity = 0.0
    hr_scenario = "no_HR"

    if dup_length >= 2 * hr_min_consec:
        inward = min(dup_length // 2, hr_max_inward)

        # Adaptive outward extension based on duplication size.
        # Large duplications get a larger outward window because flanking
        # repeats (rRNA operons, IS elements) often sit further from boundary.
        if dup_length >= 50000:
            adaptive_outward = max(hr_outward_ext, 5000)
        elif dup_length >= 10000:
            adaptive_outward = max(hr_outward_ext, 2000)
        else:
            adaptive_outward = hr_outward_ext

        # Each scenario uses strictly disjoint window halves to ensure
        # mechanistically distinct geometries: s3 outside-outside, s1
        # outside-inside, s2 inside-outside, s4 inside-inside. Module 2
        # has already refined boundaries to the precise junction position,
        # so no inside-extension tolerance on the "outside" side is needed.

        def _check_pair(seq_a, seq_b, label):
            """Run alignment + complexity filter. Returns (len, identity) or None."""
            if len(seq_a) < hr_min_consec or len(seq_b) < hr_min_consec:
                return None
            wlen, wident, _ = _run_local_alignment(
                seq_a, seq_b,
                min_window=hr_min_consec,
                min_identity=hr_min_identity,
            )
            if wlen >= hr_min_consec and wident >= hr_min_identity:
                if hr_complexity_filter:
                    if (not _hr_sequence_complexity_ok(seq_a) and
                            not _hr_sequence_complexity_ok(seq_b)):
                        return None
                return (wlen, wident, label)
            return None

        # Scenario 3: R-Y-R (repeat strictly OUTSIDE Y on both sides) — primary check
        # Boundaries are taken at face value (Module 2 has already refined them).
        # Strictly outside windows ensure s3 is mechanistically distinct from s4.
        win_start_s3 = ref_seq[max(0, start_0 - adaptive_outward) : start_0]
        win_end_s3 = ref_seq[end_0 : min(seq_len, end_0 + adaptive_outward)]

        # Scenario 1: R OUTSIDE start, INSIDE end
        # Outside-anchor side strictly outside; inside-anchor side inward.
        win_start_s1 = ref_seq[max(0, start_0 - adaptive_outward) : start_0]
        win_end_s1 = ref_seq[max(0, end_0 - inward) : end_0]

        # Scenario 2: R INSIDE start, OUTSIDE end
        # Inside-anchor side inward; outside-anchor side strictly outside.
        win_start_s2 = ref_seq[start_0 : min(seq_len, start_0 + inward)]
        win_end_s2 = ref_seq[end_0 : min(seq_len, end_0 + adaptive_outward)]

        # Scenario 4: R INSIDE Y at both boundaries
        # Windows must be disjoint, so clamp inward to half the duplication length.
        # We additionally cap at hr_max_inward (same default as other scenarios).
        s4_inward = min(inward, max(0, dup_length // 2 - 1))
        if s4_inward >= hr_min_consec:
            win_start_s4 = ref_seq[start_0 : start_0 + s4_inward]
            win_end_s4 = ref_seq[end_0 - s4_inward : end_0]
        else:
            win_start_s4 = ""
            win_end_s4 = ""

        # Try all four scenarios; keep the longest valid match
        results = []
        for seq_a, seq_b, label in [
            (win_start_s3, win_end_s3, f"single_copy_RYR_s3_w{inward}_ext{adaptive_outward}"),
            (win_start_s1, win_end_s1, f"single_copy_RYR_s1_w{inward}_ext{adaptive_outward}"),
            (win_start_s2, win_end_s2, f"single_copy_RYR_s2_w{inward}_ext{adaptive_outward}"),
            (win_start_s4, win_end_s4, f"single_copy_RYR_s4_internal_w{s4_inward}"),
        ]:
            result = _check_pair(seq_a, seq_b, label)
            if result is not None:
                results.append(result)

        if results:
            # Pick the longest match
            best = max(results, key=lambda r: r[0])
            is_hr = True
            hr_match_len = best[0]
            hr_identity = best[1]
            hr_scenario = best[2]
    else:
        hr_scenario = "too_small_for_HR_check"

    return {
        "is_hr_signature": is_hr,
        "hr_match_len": hr_match_len,
        "hr_identity": round(hr_identity, 4),
        "hr_scenario": hr_scenario,
        "microhomology_bp": microhomology_bp,
        "microhomology_seq": microhomology_seq,
        "junction_gap_bp": 0,
        "junction_seq": "",
    }


# =============================================================================
# Junction candidate generation (for modules 2 & 3)
# =============================================================================

def generate_junction_candidate(ref_seq, start_pos, end_pos, k=EXTEND_START_K):
    """Generate a junction sequence for a tandem duplication.

    Junction = ref[end-k : end] + ref[start : start+k]
    (end of copy1 followed by start of copy2)
    """
    if end_pos - k < 0 or start_pos + k > len(ref_seq):
        return None
    if end_pos > len(ref_seq) or start_pos < 0:
        return None

    left_part = ref_seq[end_pos - k : end_pos]
    right_part = ref_seq[start_pos : start_pos + k]
    return left_part + right_part


def extend_to_unique(ref_seq, start_pos, end_pos, kmer_index,
                     min_k=EXTEND_START_K, max_k=EXTEND_MAX_K):
    """Extend a junction reference until it is unique in the reference genome.

    Biological rationale:
      We only require the FULL JUNCTION SEQUENCE (left_part + right_part)
      to be absent from the reference. Individual flanks may exist in
      repetitive regions of the genome — that's fine. What makes a junction
      detectable is that its specific combination (end of region A + start
      of region B) is a novel sequence formed only by the rearrangement.

      This is particularly important for bacteria with extensive repeats
      (e.g. Streptomyces, Mycobacterium) where the flanks of duplications
      often overlap with repetitive elements, but the junction sequence
      itself is still unique.
    """
    for k in range(min_k, max_k + 1):
        junction = generate_junction_candidate(ref_seq, start_pos, end_pos, k)
        if junction is None:
            return None

        if kmer_index.is_absent(junction):
            left_part = ref_seq[end_pos - k : end_pos]
            right_part = ref_seq[start_pos : start_pos + k]

            return {
                "sequence": junction,
                "k": k,
                "left_part": left_part,
                "right_part": right_part,
                "start": start_pos,
                "end": end_pos,
                "unique": True,
            }

    return None


def generate_all_junction_candidates(ref_seq, s_approx, e_approx,
                                      window=200, flank=200,
                                      kmer_index=None, threads=1,
                                      sequences=None,
                                      coverage_array=None,
                                      genome_median=None,
                                      max_k=EXTEND_MAX_K):
    """Generate all candidate junctions for a detected region.

    Tests every (start, end) position pair in the search range.
    Uses ProcessPoolExecutor for true multi-core parallelism.

    Args:
        coverage_array, genome_median: kept for API compatibility, ignored.
        max_k: maximum half-length of junction probe (default 75 for Illumina).
    """
    import time
    from concurrent.futures import ProcessPoolExecutor, as_completed

    search_range = window + flank
    ref_len = len(ref_seq)

    s_min = max(0, s_approx - search_range)
    s_max = min(ref_len, s_approx + search_range)
    e_min = max(0, e_approx - search_range)
    e_max = min(ref_len, e_approx + search_range)

    all_starts = list(range(s_min, s_max))
    all_ends = list(range(e_min, e_max))
    total_pairs = len(all_starts) * len(all_ends)

    logger.info(
        f"  Search range: start [{s_min:,}-{s_max:,}] ({len(all_starts):,} positions), "
        f"end [{e_min:,}-{e_max:,}] ({len(all_ends):,} positions)"
    )
    logger.info(f"  Testing all {total_pairs:,} pairs")

    if not all_starts or not all_ends:
        return []

    # Junction generation with extend_to_unique
    start_time = time.time()

    if threads <= 1 or sequences is None:
        # Single-process path
        chunk_results = []
        skipped = 0
        for e_i in all_ends:
            for s_j in all_starts:
                if e_i <= s_j or (e_i - s_j) < 100:
                    continue
                result = extend_to_unique(ref_seq, s_j, e_i, kmer_index,
                                          max_k=max_k)
                if result is None:
                    skipped += 1
                    continue
                result["s_candidate"] = s_j
                result["e_candidate"] = e_i
                result["dup_size"] = e_i - s_j
                chunk_results.append(result)
        all_results = chunk_results
        total_skipped = skipped
    else:
        # Multi-process: split filtered end positions across processes
        # Use 10× more chunks than cores for load balancing.
        # ProcessPoolExecutor assigns next chunk to whichever core
        # finishes first — slow repetitive chunks don't block fast ones.
        n_chunks = min(threads * 10, len(all_ends))
        if n_chunks == 0:
            return []

        chunks = []
        for i in range(n_chunks):
            chunk_start = i * len(all_ends) // n_chunks
            chunk_end = (i + 1) * len(all_ends) // n_chunks
            chunk_positions = all_ends[chunk_start:chunk_end]
            chunks.append((sequences, ref_seq, chunk_positions,
                          all_starts, max_k))

        logger.info(
            f"  Split into {len(chunks)} chunks across {threads} processes"
        )

        all_results = []
        total_skipped = 0
        completed = 0

        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {
                executor.submit(_process_end_position_chunk_filtered, chunk): i
                for i, chunk in enumerate(chunks)
            }
            for future in as_completed(futures):
                chunk_results, chunk_skipped = future.result()
                all_results.extend(chunk_results)
                total_skipped += chunk_skipped
                completed += 1

                # Log progress every ~5%
                log_interval = max(1, len(chunks) // 20)
                if completed % log_interval == 0 or completed == len(chunks):
                    elapsed = time.time() - start_time
                    pct = 100 * completed / len(chunks)
                    logger.info(
                        f"    Progress: {completed}/{len(chunks)} chunks ({pct:.0f}%), "
                        f"{len(all_results)} candidates so far, {elapsed:.1f}s"
                    )

    # Deduplicate by junction sequence
    seen_sequences = set()
    unique_junctions = []
    for result in all_results:
        seq = result["sequence"]
        if seq not in seen_sequences:
            seen_sequences.add(seq)
            unique_junctions.append(result)

    elapsed = time.time() - start_time
    logger.info(
        f"  Found {len(unique_junctions)} unique junction candidates "
        f"(skipped {total_skipped} non-unique, {elapsed:.1f}s)"
    )
    return unique_junctions


def _process_end_position_chunk_filtered(args):
    """Worker function for parallel junction generation with filtered positions.

    Takes filtered start positions (list) instead of a range.
    """
    sequences, ref_seq, e_positions, s_positions, max_k = args

    # Each worker builds its own KmerIndex
    worker_index = KmerIndex(sequences)

    results = []
    skipped = 0

    for e_i in e_positions:
        for s_j in s_positions:
            if e_i <= s_j or (e_i - s_j) < 100:
                continue

            result = extend_to_unique(ref_seq, s_j, e_i, worker_index, max_k=max_k)
            if result is None:
                skipped += 1
                continue

            result["s_candidate"] = s_j
            result["e_candidate"] = e_i
            result["dup_size"] = e_i - s_j
            results.append(result)

    return results, skipped


def build_junction_reference(ref_seq, start, end, kmer_index,
                              provided_sequence=None):
    """Build a junction reference for a known duplication (module 3)."""
    if provided_sequence:
        if kmer_index.is_absent(provided_sequence):
            return {
                "sequence": provided_sequence,
                "k": len(provided_sequence) // 2,
                "start": start,
                "end": end,
                "unique": True,
                "source": "user_provided",
            }
        else:
            logger.warning(
                f"  User-provided junction at ({start}, {end}) "
                f"not unique. Auto-generating."
            )

    result = extend_to_unique(ref_seq, start, end, kmer_index)
    if result:
        result["source"] = "auto_generated"
    return result
