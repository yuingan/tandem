"""Find duplication junctions from the reads, in one pass over the FASTQ.

A read crossing a tandem junction anchors to two diagonals on the reference:
the second sits behind the first by the duplication size. No coordinates
needed from the user.

Support is counted the same way on both alleles - reads crossing the junction
against reads running straight through the same boundary, same anchor length,
same tolerance - so the losses cancel:

    dosage = junction_reads / reference_reads ~ f(n-1)

A ratio, not a proportion: J/(J+R) is f/(1+f) and caps at 0.5.
"""

import gzip
import logging

import numpy as np

from .kmer_index import KmerIndex, _encode_all, reverse_complement

logger = logging.getLogger("tandem")

ANCHOR_K = 20        # k-mer length used to place a read on the reference
STRIDE = 8           # spacing between anchors along a read
MIN_ANCHOR = 25      # bases required either side of a junction
MIN_SUPPORT = 2      # reads needed before a junction is reported
CHUNK_READS = 100_000


def _open(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_chunks(paths, chunk=CHUNK_READS):
    """Yield lists of read sequences from one or more FASTQ files."""
    batch = []
    for path in paths:
        if path is None:
            continue
        with _open(path) as f:
            for i, line in enumerate(f):
                if i % 4 != 1:
                    continue
                batch.append(line.strip().upper())
                if len(batch) >= chunk:
                    yield batch
                    batch = []
    if batch:
        yield batch


def _anchor_reads(reads, index, k, stride):
    """Returns (read_id, offset in read, diagonal) per matched anchor."""
    spacer = "N" * k
    concat = spacer.join(reads)
    starts = np.cumsum([0] + [len(r) + k for r in reads[:-1]])

    codes, offsets = _encode_all(concat, k, stride=stride)
    if len(codes) == 0:
        return (np.empty(0, np.int64),) * 3

    lo = np.searchsorted(index._codes, codes, side="left")
    hi = np.searchsorted(index._codes, codes, side="right")
    # a repeated anchor tells us nothing about where the read sits
    unique = (hi - lo) == 1
    if not unique.any():
        return (np.empty(0, np.int64),) * 3

    offsets = offsets[unique]
    ref_pos = index._offsets[lo[unique]]

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
         min_anchor=MIN_ANCHOR, min_support=MIN_SUPPORT, index=None):
    """Junctions with support counts. Coordinates are 1-based inclusive."""
    if index is None:
        index = KmerIndex({"ref": ref_seq}, min_k=k)

    votes = {}
    linear_diag = []
    linear_len = []
    depth_delta = np.zeros(len(ref_seq) + 2, dtype=np.int32)

    n_reads = 0
    for reads in read_chunks(read_paths):
        n_reads += len(reads)
        for orientation in (0, 1):
            batch = reads if orientation == 0 else [reverse_complement(r) for r in reads]
            read_id, in_read, diag = _anchor_reads(batch, index, k, stride)
            if len(read_id) == 0:
                continue

            lengths = np.fromiter((len(r) for r in batch), dtype=np.int64,
                                  count=len(batch))
            edges = np.r_[0, np.nonzero(np.diff(read_id))[0] + 1]
            lo_diag = np.minimum.reduceat(diag, edges)
            hi_diag = np.maximum.reduceat(diag, edges)

            # one diagonal = ancestral allele + depth signal; handled in bulk
            straight = lo_diag == hi_diag
            d = lo_diag[straight]
            n = lengths[read_id[edges[straight]]]
            linear_diag.append(d)
            linear_len.append(n)

            starts = np.clip(d, 0, len(ref_seq))
            ends = np.clip(d + n, 0, len(ref_seq))
            np.add.at(depth_delta, starts, 1)
            np.add.at(depth_delta, ends, -1)

            for i in np.nonzero(~straight)[0]:
                lo = edges[i]
                hi = edges[i + 1] if i + 1 < len(edges) else len(read_id)
                _vote(votes, batch[read_id[lo]], ref_seq,
                      diag[lo:hi], in_read[lo:hi], k, min_anchor)

    logger.info(f"  Scanned {n_reads:,} reads, {len(votes)} candidate junction(s)")

    diags = np.concatenate(linear_diag) if linear_diag else np.empty(0, np.int64)
    lens = np.concatenate(linear_len) if linear_len else np.empty(0, np.int64)
    depth = np.cumsum(depth_delta[:len(ref_seq)])

    results = []
    for (start0, end0), support in sorted(votes.items(), key=lambda kv: -kv[1]):
        if support < min_support:
            continue
        ref_count = _count_reference_reads(diags, lens, end0 + 1, min_anchor)
        microhomology = _homology(start0, end0, ref_seq)
        results.append({
            "start": start0 + 1,
            "end": end0 + 1,
            "size": end0 - start0 + 1,
            "microhomology": microhomology,
            "junction_reads": int(support),
            "reference_reads": int(ref_count),
            "dosage": round(support / ref_count, 4) if ref_count else None,
            "depth_dosage": _depth_dosage(depth, start0, end0),
        })
    return results


def _vote(votes, read, ref_seq, diags, in_read, k, min_anchor):
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
    if split < min_anchor or len(read) - split < min_anchor:
        return

    end0 = split + uniq[0] - 1     # last base of the unit, 0-based
    start0 = split + uniq[1]       # first base of the unit, 0-based
    if start0 >= end0:
        return

    start0, end0, _ = left_align(start0, end0, ref_seq)
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


def _depth_dosage(depth, start0, end0):
    """Same quantity as the junction ratio, from coverage instead.

    Depth inside a duplicated interval is 1 + f(n-1) times the depth outside,
    so the two estimates should agree.
    """
    span = end0 - start0 + 1
    margin = 200                    # step off the junction before sampling
    inside = depth[start0 + margin:end0 + 1 - margin]
    left = depth[max(0, start0 - span - margin):max(0, start0 - margin)]
    right = depth[end0 + 1 + margin:end0 + 1 + margin + span]
    outside = np.concatenate([left, right])
    if len(inside) < margin or len(outside) < margin:
        return None
    # median, not mean: repeats in either window drop to zero depth
    base = float(np.median(outside))
    if base == 0:
        return None
    return round(float(np.median(inside)) / base - 1.0, 4)


def _count_reference_reads(diags, lengths, boundary_0, min_anchor):
    """The denominator: ancestral-boundary reads, counted like the numerator."""
    if not len(diags):
        return 0
    left = boundary_0 - diags
    return int(((left >= min_anchor) & (lengths - left >= min_anchor)).sum())
