"""Sorted-array k-mer index over the reference.

Each k-mer is packed into an int64 (2 bits per base, so k <= 31) and the
array sorted once. Lookup is a binary search plus a string compare at the
candidate positions. Replaces a str.find scan over the whole genome per query.
"""

import logging

import numpy as np

logger = logging.getLogger("tandem")

INDEX_K = 20            # must be <= 31 to fit two bits per base in an int64
SPACER = "N" * 100      # keeps k-mers from spanning two sequences

_CODE = np.full(256, -1, dtype=np.int64)
for _i, _b in enumerate("ACGT"):
    _CODE[ord(_b)] = _i
    _CODE[ord(_b.lower())] = _i

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def _rc_codes(codes, k):
    """Reverse complement of packed k-mers (A=0 C=1 G=2 T=3, last base lowest)."""
    rc = np.zeros_like(codes)
    x = codes.copy()
    for _ in range(k):
        rc = (rc << 2) | (3 - (x & 3))
        x >>= 2
    return rc


def _encode_all(seq, k, stride=1):
    """Pack every stride-th k-mer into an int64. K-mers with N are dropped."""
    values = _CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]
    n = len(values) - k + 1
    if n <= 0:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)

    offsets = np.arange(0, n, stride, dtype=np.int64)
    codes = np.zeros(len(offsets), dtype=np.int64)
    valid = np.ones(len(offsets), dtype=bool)
    for i in range(k):
        window = values[offsets + i]
        codes = (codes << 2) | np.maximum(window, 0)
        valid &= window >= 0

    keep = np.nonzero(valid)[0]
    return codes[keep], offsets[keep]


class KmerIndex:
    """count_occurrences / is_unique / is_absent, plus positions() for the scanner."""

    def __init__(self, sequences, min_k=INDEX_K):
        self.k = min(min_k, 31)
        self.sequences = sequences

        parts = []
        self._seq_offsets = {}
        offset = 0
        for seq_id, seq in sequences.items():
            self._seq_offsets[seq_id] = offset
            parts.append(seq.upper())
            parts.append(SPACER)
            offset += len(seq) + len(SPACER)
        self._concat = "".join(parts)

        codes, offsets = _encode_all(self._concat, self.k)
        order = np.argsort(codes, kind="stable")
        self._codes = codes[order]
        self._offsets = offsets[order]
        self._anchors = None

        logger.debug(
            f"KmerIndex: {len(sequences)} sequence(s), "
            f"{len(self._concat):,} bp, {len(self._codes):,} k-mers at k={self.k}"
        )

    def anchor_table(self):
        """K-mers that can place a read, sorted, with their positions.

        An anchor occurs once in the genome and not at all on the other
        strand. Working this out once here leaves the scanner one binary
        search per read k-mer instead of four.
        """
        if self._anchors is None:
            codes, offsets = self._codes, self._offsets
            if len(codes) == 0:
                self._anchors = (codes, offsets)
                return self._anchors
            differs = codes[1:] != codes[:-1]
            once = np.r_[True, differs] & np.r_[differs, True]
            c, p = codes[once], offsets[once]
            rc = _rc_codes(c, self.k)
            j = np.minimum(np.searchsorted(codes, rc), len(codes) - 1)
            keep = codes[j] != rc            # palindromes fall out here too
            self._anchors = (c[keep], p[keep])
        return self._anchors

    def anchors_at(self, codes):
        """Index into anchor_table() for each code, or -1 where it is no anchor."""
        table, _ = self.anchor_table()
        if len(table) == 0 or len(codes) == 0:
            return np.full(len(codes), -1, dtype=np.int64)
        j = np.minimum(np.searchsorted(table, codes), len(table) - 1)
        return np.where(table[j] == codes, j, -1)

    def _code(self, kmer):
        values = _CODE[np.frombuffer(kmer.encode(), dtype=np.uint8)]
        if (values < 0).any():
            return None
        code = 0
        for v in values:
            code = (code << 2) | int(v)
        return code

    def _seed_positions(self, kmer):
        code = self._code(kmer)
        if code is None:
            return np.empty(0, dtype=np.int64)
        lo = np.searchsorted(self._codes, code, side="left")
        hi = np.searchsorted(self._codes, code, side="right")
        return self._offsets[lo:hi]

    def positions(self, query):
        """Offsets in the concatenated reference where query occurs, forward only."""
        q = query.upper()
        if len(q) < self.k:
            return self._scan(q)

        seed = q[:self.k]
        hits = []
        for pos in self._seed_positions(seed):
            pos = int(pos)
            if self._concat[pos:pos + len(q)] == q:
                hits.append(pos)
        return hits

    def _scan(self, q):
        """For queries shorter than k."""
        hits = []
        start = 0
        while True:
            pos = self._concat.find(q, start)
            if pos == -1:
                return hits
            hits.append(pos)
            start = pos + 1

    def count_occurrences(self, query, max_count=None):
        """Exact occurrences on both strands."""
        q = query.upper()
        count = len(self.positions(q))
        if max_count is not None and count >= max_count:
            return count

        rc = reverse_complement(q)
        if rc != q:
            count += len(self.positions(rc))
        return count

    def is_unique(self, query):
        return self.count_occurrences(query, max_count=2) <= 1

    def is_absent(self, query):
        return self.count_occurrences(query, max_count=1) == 0

    def locate(self, offset):
        """Concatenated offset -> (seq_id, position)."""
        best = None
        for seq_id, start in self._seq_offsets.items():
            if offset >= start and (best is None or start > best[1]):
                best = (seq_id, start)
        if best is None:
            return None, offset
        return best[0], offset - best[1]
