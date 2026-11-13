"""Find duplication junctions from the reads, in one pass over the FASTQ.

A read crossing a tandem junction anchors to two diagonals on the reference:
the second sits behind the first by the duplication size. No coordinates
needed from the user.

Support is counted the same way on both alleles - reads crossing the junction
against reads running straight through the same boundary, same anchor length,
same tolerance - so the losses cancel:

    dosage = junction_reads / reference_reads ~ f(n-1)

A ratio, not a proportion: J/(J+R) is f/(1+f) and caps at 0.5.

Both kinds of read count only with a unique anchor wholly on each side of
the boundary they cross. Otherwise a junction read whose far side sits in a
repeat keeps anchors on one side only and passes for a reference read.

Next to a repeat the two counts still lose different shares of their reads,
so each is divided by the number of read placements the scan could have
counted at that spot, worked out from the reference (_countable). The
uncorrected ratio is kept as dosage_raw.

dosage_low and dosage_high bound a 95% interval: for the junction ratio from
the two counts, for depth from how much the same depth statistic moves where
nothing is duplicated.
"""

import gzip
import logging
import math
import multiprocessing
import time
from collections import deque

import numpy as np

from .kmer_index import KmerIndex, _encode_all, reverse_complement

logger = logging.getLogger("tandem")

ANCHOR_K = 20        # k-mer length used to place a read on the reference
STRIDE = 8           # spacing between anchors along a read
MIN_ANCHOR = 25      # bases required either side of a junction
MIN_SUPPORT = 2      # reads needed before a junction is reported
MIN_SIZE = 200       # smaller "duplications" are mostly indels; same floor as module 1
ORIGIN_SLOP = 1000   # see _vote
CHUNK_READS = 100_000


def _open(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


class FastqError(ValueError):
    pass


def read_chunks(paths, chunk=CHUNK_READS):
    """Yield lists of read sequences from one or more FASTQ files.

    Every record is checked on the way past: a header starting with @, a
    separator starting with +, as many quality values as bases. One broken
    record would otherwise put every later line out of step, and the rest of
    the file would be lost without a word.
    """
    batch = []
    for path in paths:
        if path is None:
            continue
        with _open(path) as f:
            part, seq = 0, ""
            for n, line in enumerate(f, 1):
                if part == 0:
                    if not line.strip():
                        continue                # blank lines between records
                    if not line.startswith("@"):
                        raise FastqError(f"{path}, line {n}: expected a read header "
                                         f"starting with @, found {line[:40].rstrip()!r}")
                elif part == 1:
                    seq = line.strip()
                elif part == 2:
                    if not line.startswith("+"):
                        raise FastqError(f"{path}, line {n}: expected the + line of a "
                                         f"read, found {line[:40].rstrip()!r}")
                else:
                    if len(line.strip()) != len(seq):
                        raise FastqError(f"{path}, line {n}: {len(seq)} bases but "
                                         f"{len(line.strip())} quality values")
                    batch.append(seq.upper())
                    if len(batch) >= chunk:
                        yield batch
                        batch = []
                part = (part + 1) % 4
            if part:
                raise FastqError(f"{path}: the last read is incomplete")
    if batch:
        yield batch


def _anchor_reads(reads, index, k, stride, contig_start, contig_end):
    """Returns (read_id, offset in read, diagonal) per matched anchor.

    Only anchors on [contig_start, contig_end) of the concatenated index are
    kept, and positions come back in that contig's own coordinates.
    """
    # Pad each read so the next one starts on a multiple of the stride. Every
    # read then has its anchors at offsets 0, stride, 2*stride..., the one
    # pattern _countable assumes, whatever the read lengths are.
    pads = [k + (-(len(r) + k)) % stride for r in reads]
    concat = "".join(r + "N" * pad for r, pad in zip(reads, pads))
    starts = np.cumsum([0] + [len(r) + pad for r, pad in zip(reads[:-1], pads[:-1])])

    codes, offsets = _encode_all(concat, k, stride=stride)
    if len(codes) == 0:
        return (np.empty(0, np.int64),) * 3

    # only k-mers found once in the genome, and not on the other strand, can
    # place a read; the index has that list ready
    hit = index.anchors_at(codes)
    keep = hit >= 0
    if not keep.any():
        return (np.empty(0, np.int64),) * 3

    _, anchor_pos = index.anchor_table()
    offsets = offsets[keep]
    ref_pos = anchor_pos[hit[keep]]
    on_contig = (ref_pos >= contig_start) & (ref_pos < contig_end)
    offsets = offsets[on_contig]
    ref_pos = ref_pos[on_contig] - contig_start
    if not len(offsets):
        return (np.empty(0, np.int64),) * 3

    read_id = np.searchsorted(starts, offsets, side="right") - 1
    in_read = offsets - starts[read_id]
    return read_id, in_read, ref_pos - in_read


SPLIT_WINDOW = 30    # bases either side of the junction used to place it
MAX_SPLIT_MISMATCH = 3


def _exact_junction(read, ref, diag_left, diag_right, lo, hi):
    """Offset in the read where the left diagonal gives way to the right one.

    Only bases near the junction are compared, so an error elsewhere in the
    read cannot decide the answer.
    """
    n = len(read)
    best, best_cost = None, None

    for split in range(max(lo, 0), min(hi, n) + 1):
        cost = 0
        for x in range(max(0, split - SPLIT_WINDOW), split):
            p = x + diag_left
            if p < 0 or p >= len(ref) or ref[p] != read[x]:
                cost += 1
        for x in range(split, min(n, split + SPLIT_WINDOW)):
            p = x + diag_right
            if p < 0 or p >= len(ref) or ref[p] != read[x]:
                cost += 1
        if best_cost is None or cost < best_cost:
            best, best_cost = split, cost

    if best_cost is None or best_cost > MAX_SPLIT_MISMATCH:
        return None
    return best


def scan(ref_seq, read_paths, k=ANCHOR_K, stride=STRIDE,
         min_anchor=MIN_ANCHOR, min_support=MIN_SUPPORT, min_size=MIN_SIZE,
         index=None, contig_start=0, threads=1, targets=None, slop=20):
    """Junctions with support counts. Coordinates are 1-based inclusive.

    With a shared index over several contigs, contig_start is where ref_seq
    begins in that index (KmerIndex._seq_offsets). With threads > 1 blocks of
    reads go to worker processes; results are merged in file order, so the
    output does not depend on the thread count.

    targets: (start, end) of known duplicated units. Instead of the junctions
    found, each target is measured, in order, whatever its support: at the
    best-supported junction within slop bases of both ends, or where it is
    if no read crosses it.
    """
    if index is None:
        index = KmerIndex({"ref": ref_seq}, min_k=k)
        contig_start = 0
    k = index.k
    index.anchor_table()          # build once here, before any worker forks
    _SHARED.update(index=index, ref_seq=ref_seq, k=k, stride=stride,
                   contig_start=contig_start, contig_end=contig_start + len(ref_seq),
                   min_anchor=min_anchor, min_size=min_size)

    votes = {}
    parts = {"diag": [], "len": [], "first": [], "last": []}
    depth_delta = np.zeros(len(ref_seq) + 2, dtype=np.int32)
    n_reads = 0
    next_report = 1_000_000
    started = time.time()

    def merge(result):
        nonlocal n_reads, next_report
        n, chunk_votes, arrays = result
        n_reads += n
        for key, v in chunk_votes.items():
            votes[key] = votes.get(key, 0) + v
        for key in parts:
            parts[key].append(arrays[key])
        np.add.at(depth_delta, np.clip(arrays["diag"], 0, len(ref_seq)), 1)
        np.add.at(depth_delta, np.clip(arrays["diag"] + arrays["len"], 0, len(ref_seq)), -1)
        if n_reads >= next_report:
            logger.info(f"  {n_reads:,} reads scanned ({time.time() - started:.0f}s)")
            next_report += 1_000_000

    if threads > 1:
        # fork, so workers share the index rather than each building one
        with multiprocessing.get_context("fork").Pool(threads) as pool:
            pending = deque()
            for reads in read_chunks(read_paths):
                pending.append(pool.apply_async(_scan_chunk, (reads,)))
                if len(pending) >= 2 * threads:     # keep memory bounded
                    merge(pending.popleft().get())
            while pending:
                merge(pending.popleft().get())
    else:
        for reads in read_chunks(read_paths):
            merge(_scan_chunk(reads))

    logger.info(f"  Scanned {n_reads:,} reads, {len(votes)} candidate junction(s)")

    empty = np.empty(0, np.int64)
    diags, lens, firsts, lasts = (np.concatenate(parts[key]) if parts[key] else empty
                                  for key in ("diag", "len", "first", "last"))
    depth = np.cumsum(depth_delta[:len(ref_seq)])
    read_len = int(np.median(lens)) if len(lens) else 150
    logger.info(f"  Median read length {read_len} bp")

    def measure(start0, end0, support):
        ref_count = _count_reference_reads(diags, lens, firsts, lasts,
                                           end0 + 1, min_anchor, k)
        mh = _homology(start0, end0, ref_seq)
        uniq = _uniqueness(ref_seq, index, k, [(end0 - read_len, end0 + read_len),
                                               (start0 - read_len, start0 + read_len)])
        can_j = _countable(uniq, end0, start0, read_len, min_anchor, k, stride, mh)
        can_r = _countable(uniq, end0, end0 + 1, read_len, min_anchor, k, stride)
        most = read_len - 2 * min_anchor + 1     # placements with no repeats around

        raw = support / ref_count if ref_count else None
        fixed = raw * can_r / can_j if raw is not None and can_j else None
        from_depth, noise = _depth_dosage(
            depth, start0, end0, unique=lambda a, b: _unique_mask(ref_seq, index, k, a, b))
        # When no read can be placed on both sides of the junction (one end of
        # the unit inside a repeat longer than a read) or of the reference
        # boundary, the reads can't give a ratio; coverage still can.
        dosage, source, low, high = fixed, "junction", None, None
        if fixed is not None:
            low, high = (x * can_r / can_j for x in _ratio_interval(support, ref_count))
        elif from_depth is not None:
            dosage, source = from_depth, "depth"
            if noise is not None:
                low, high = from_depth - 2 * noise, from_depth + 2 * noise
        else:
            source = ""
        return {
            "start": start0 + 1,
            "end": end0 + 1,
            "size": end0 - start0 + 1,
            "microhomology": mh,
            "junction_reads": int(support),
            "reference_reads": int(ref_count),
            "dosage": _round(dosage),
            "dosage_from": source,
            "dosage_low": _round(low),
            "dosage_high": _round(high),
            "dosage_raw": _round(raw),
            "junction_visible": round(can_j / most, 3),
            "reference_visible": round(can_r / most, 3),
            "depth_dosage": from_depth,
            "depth_noise": noise,
        }

    if targets is not None:
        measured = []
        for start, end in targets:
            start0, end0, _ = left_align(start - 1, end - 1, ref_seq)
            near = [(n, -abs(s - start0) - abs(e - end0), s, e)
                    for (s, e), n in votes.items()
                    if abs(s - start0) <= slop and abs(e - end0) <= slop]
            support = 0
            if near:
                support, _, start0, end0 = max(near)
            measured.append(measure(start0, end0, support))
        return measured

    return [measure(start0, end0, support)
            for (start0, end0), support in sorted(votes.items(), key=lambda kv: -kv[1])
            if support >= min_support]


_SHARED = {}    # what _scan_chunk needs; set by scan() before workers fork


def _scan_chunk(reads):
    """Anchor one block of reads (both strands) and vote for junctions.

    Returns the read count, the votes, and the geometry of the reads that sit
    on a single diagonal, which scan() needs for the reference count and depth.
    """
    sh = _SHARED
    k, stride, ref_seq = sh["k"], sh["stride"], sh["ref_seq"]
    votes = {}
    out = {"diag": [], "len": [], "first": [], "last": []}

    for orientation in (0, 1):
        batch = reads if orientation == 0 else [reverse_complement(r) for r in reads]
        read_id, in_read, diag = _anchor_reads(batch, sh["index"], k, stride,
                                               sh["contig_start"], sh["contig_end"])
        if len(read_id) == 0:
            continue

        lengths = np.fromiter((len(r) for r in batch), dtype=np.int64, count=len(batch))
        edges = np.r_[0, np.nonzero(np.diff(read_id))[0] + 1]
        lo_diag = np.minimum.reduceat(diag, edges)
        hi_diag = np.maximum.reduceat(diag, edges)
        first = np.minimum.reduceat(in_read, edges)
        last = np.maximum.reduceat(in_read, edges)

        # one diagonal = ancestral allele + depth signal; handled in bulk
        straight = lo_diag == hi_diag
        out["diag"].append(lo_diag[straight])
        out["len"].append(lengths[read_id[edges[straight]]])
        out["first"].append(first[straight])
        out["last"].append(last[straight])

        for i in np.nonzero(~straight)[0]:
            lo = edges[i]
            hi = edges[i + 1] if i + 1 < len(edges) else len(read_id)
            _vote(votes, batch[read_id[lo]], ref_seq, diag[lo:hi], in_read[lo:hi],
                  k, sh["min_anchor"], sh["min_size"])

    empty = np.empty(0, np.int64)
    arrays = {key: (np.concatenate(v) if v else empty) for key, v in out.items()}
    return len(reads), votes, arrays


def _vote(votes, read, ref_seq, diags, in_read, k, min_anchor, min_size):
    """Record where one read puts the junction."""
    # second diagonal behind the first = tandem; anything else is some other
    # rearrangement and not ours to call
    uniq = list(dict.fromkeys(int(v) for v in diags))
    if len(uniq) != 2 or uniq[1] >= uniq[0]:
        return

    on_left = in_read[diags == uniq[0]]
    on_right = in_read[diags == uniq[1]]
    split = _exact_junction(read, ref_seq, uniq[0], uniq[1],
                            int(on_left.max()) + k, int(on_right.min()))
    if split is None:
        return

    end0 = split + uniq[0] - 1     # last base of the unit, 0-based
    start0 = split + uniq[1]       # first base of the unit, 0-based
    if start0 >= end0:
        return
    start0, end0, shift = left_align(start0, end0, ref_seq)
    # enough bases on each side, measured where the junction is reported
    s = split - shift
    if s < min_anchor or len(read) - s < min_anchor:
        return
    # A read across the origin of a circular contig looks exactly like a
    # duplication of the whole contig. The slop covers assemblies whose ends
    # overlap by a few hundred bases.
    if start0 <= ORIGIN_SLOP and end0 >= len(ref_seq) - 1 - ORIGIN_SLOP:
        return
    if end0 - start0 + 1 < min_size:
        return
    votes[(start0, end0)] = votes.get((start0, end0), 0) + 1


def left_align(start0, end0, ref):
    """Shift a junction to the leftmost position giving the same sequence.

    Without this, reads settle on different positions within the microhomology
    and the support for one junction splits across several entries.
    """
    shift = 0
    while start0 > 0 and ref[end0] == ref[start0 - 1]:
        start0 -= 1
        end0 -= 1
        shift += 1
    return start0, end0, shift


def _homology(start0, end0, ref):
    """Microhomology: how far a left-aligned junction can slide right."""
    n = 0
    while end0 + 1 + n < len(ref) and ref[start0 + n] == ref[end0 + 1 + n]:
        n += 1
    return n


def _round(x, digits=4):
    if x is None:
        return None
    return round(x, digits) if math.isfinite(x) else x


def _ratio_interval(j, r, z=1.96):
    """95% interval for j/r, from a Wilson interval on j's share of j + r."""
    n = j + r
    p = j / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    lo = 0.0 if j == 0 else max(0.0, centre - half)
    hi = min(1.0, centre + half)
    return lo / (1 - lo), (hi / (1 - hi) if hi < 1 else math.inf)


def _centre(values):
    """Mean of the middle half: robust like a median, but not stuck on whole
    numbers, which would make neighbouring windows tie and the noise read 0."""
    lo, hi = np.percentile(values, [25, 75])
    middle = values[(values >= lo) & (values <= hi)]
    return float(middle.mean()) if len(middle) else float(np.median(values))


def _depth_dosage(depth, start0, end0, unique=None, windows=30, longest=50_000):
    """Same quantity as the junction ratio, from coverage instead, and its noise.

    Depth inside a duplicated interval is 1 + f(n-1) times the depth outside,
    so the two estimates should agree. Only positions whose k-mer is unique in
    the genome are used (unique(a, b) gives the mask): few reads can be placed
    in a repeat, and a repeat in the flanks would otherwise lower the baseline
    and inflate the estimate.

    The noise is how far the same statistic strays where nothing is
    duplicated: unit-sized windows are laid end to end beyond the baseline, up
    to `windows` on each side, and each is set against its two neighbours as
    the unit is against its flanks. On simulated coverage it matches the
    spread of the estimate itself from 1.5 to 45 kb.
    """
    def take(a, b):
        a, b = max(a, 0), min(b, len(depth))
        values = depth[a:b]
        if unique is None or b <= a:
            return values
        keep = np.zeros(b - a, dtype=bool)
        mask = unique(a, b)
        keep[:len(mask)] = mask[:b - a]
        return values[keep]

    span = end0 - start0 + 1
    margin = 200                    # step off the junction before sampling
    inside = take(start0 + margin, end0 + 1 - margin)
    left_from = max(0, start0 - span - margin)
    right_to = end0 + 1 + margin + span
    outside = np.concatenate([take(left_from, start0 - margin),
                              take(end0 + 1 + margin, right_to)])
    if len(inside) < margin or len(outside) < margin:
        return None, None
    base = _centre(outside)
    if base == 0:
        return None, None
    dosage = round(_centre(inside) / base - 1.0, 4)

    w = min(span - 2 * margin, longest)
    null = []
    for starts in ([left_from - (k + 1) * w for k in range(windows)],
                   [right_to + k * w for k in range(windows)]):
        tiles = [take(a, a + w) for a in starts if 0 <= a and a + w <= len(depth)]
        tiles = [t for t in tiles if len(t) >= margin]
        for prev, tile, nxt in zip(tiles, tiles[1:], tiles[2:]):
            flank = _centre(np.concatenate([prev, nxt]))
            if flank > 0:
                null.append(_centre(tile) / flank - 1.0)
    if len(null) < 5:
        return dosage, None
    # robust spread: another duplication among the windows shouldn't count
    mad = float(np.median(np.abs(np.array(null) - np.median(null))))
    return dosage, round(1.4826 * mad, 4)


def _unique_mask(ref_seq, index, k, a, b):
    """For positions a..b-1 (a >= 0): does the k-mer starting there occur once?"""
    b = min(b, len(ref_seq) - k + 1)
    mask = np.zeros(max(b - a, 0), dtype=bool)
    if b > a:
        codes, offs = _encode_all(ref_seq[a:b + k - 1], k)
        mask[offs] = index.anchors_at(codes) >= 0
    return mask


def _uniqueness(ref_seq, index, k, windows):
    """uniq(p): does the k-mer starting at p occur once in the genome, either strand?"""
    masks = []
    for a, b in windows:
        a = max(a, 0)
        masks.append((a, _unique_mask(ref_seq, index, k, a, b)))

    def uniq(p):
        for a, mask in masks:
            if a <= p < a + len(mask):
                return bool(mask[p - a])
        return False
    return uniq


def _countable(uniq, left_end, right_start, read_len, min_anchor, k, stride, mh=0):
    """Read placements across a boundary the scan would count.

    The read takes s bases ending at left_end, then the rest from right_start
    (left_end + 1 for plain reference). Its anchors sit at offsets 0, stride,
    2*stride... (see _anchor_reads), and it counts with a unique anchor wholly
    on each side, the rule both _vote and _count_reference_reads apply.

    At a junction with mh bases of microhomology, s is measured at the
    left-aligned position, and the first mh bases after it read the same on
    both sides, so an anchor ending inside them still places the left side.
    """
    total = 0
    for s in range(min_anchor, read_len - min_anchor + 1):
        left = right = False
        for o in range(0, read_len - k + 1, stride):
            if o + k <= s + mh:
                left = left or uniq(left_end - s + 1 + o)
            elif o >= s:
                right = right or uniq(right_start + o - s)
        total += left and right
    return total


def _count_reference_reads(diags, lengths, firsts, lasts, boundary_0, min_anchor, k):
    """The denominator: reads straight through the ancestral boundary.

    Same rule as a junction read: enough bases on each side, and a unique
    anchor wholly on each side.
    """
    if not len(diags):
        return 0
    left = boundary_0 - diags            # bases of the read before the boundary
    ok = ((left >= min_anchor) & (lengths - left >= min_anchor)
          & (firsts + k <= left) & (lasts >= left))
    return int(ok.sum())
