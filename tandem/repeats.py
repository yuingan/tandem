"""Count how many places in the genome carry each candidate copy.

A copy that also appears at a dozen other loci is an IS element or an rRNA
operon, not a one-off duplication. The self-alignment already contains every
one of those loci, so this is just bookkeeping over it.
"""

import logging
from bisect import bisect_right

logger = logging.getLogger("tandem")

MIN_OVERLAP = 0.5   # fraction of the query a hit must cover to count


def _span(start, end):
    return (min(start, end), max(start, end))


def _overlap(a, b):
    return min(a[1], b[1]) - max(a[0], b[0]) + 1


def _same_locus(a, b, reciprocal=0.9):
    """Same place, not just touching (tandem copies share a base at the junction)."""
    shared = _overlap(a, b)
    if shared <= 0:
        return False
    return (shared >= (a[1] - a[0] + 1) * reciprocal
            and shared >= (b[1] - b[0] + 1) * reciprocal)


class LocusIndex:
    """Lookup from an interval to the other loci its sequence aligns to."""

    def __init__(self, alignments):
        raw = {}
        for aln in alignments:
            ref = _span(aln["ref_start"], aln["ref_end"])
            qry = _span(aln["qry_start"], aln["qry_end"])
            # skip the self-diagonal, it overlaps everything
            if aln["ref_id"] == aln["qry_id"] and _same_locus(ref, qry):
                continue
            raw.setdefault(aln["ref_id"], []).append((ref, aln["qry_id"], qry))
            raw.setdefault(aln["qry_id"], []).append((qry, aln["ref_id"], ref))

        self.index = {}
        for seq_id, entries in raw.items():
            entries.sort(key=lambda e: e[0][0])
            starts = [e[0][0] for e in entries]
            longest = max(e[0][1] - e[0][0] for e in entries)
            self.index[seq_id] = (entries, starts, longest)

    def partners(self, seq_id, start, end):
        """Loci aligning over at least MIN_OVERLAP of [start, end]."""
        if seq_id not in self.index:
            return []

        entries, starts, longest = self.index[seq_id]
        query = (start, end)
        needed = (end - start + 1) * MIN_OVERLAP

        lo = bisect_right(starts, start - longest - 1)
        hi = bisect_right(starts, end)

        hits = []
        for i in range(lo, hi):
            interval, partner_seq, partner_interval = entries[i]
            if _overlap(query, interval) >= needed:
                hits.append((partner_seq, partner_interval))
        return hits


def _merge(intervals):
    merged = []
    for interval in sorted(intervals):
        if merged and _same_locus(merged[-1], interval, reciprocal=MIN_OVERLAP):
            last = merged[-1]
            merged[-1] = (min(last[0], interval[0]), max(last[1], interval[1]))
            continue
        merged.append(interval)
    return merged


def count_loci(index, seq_id, start, end):
    """Loci carrying this sequence, query included: unique = 1, tandem pair = 2."""
    by_seq = {}
    for partner_seq, interval in index.partners(seq_id, start, end):
        by_seq.setdefault(partner_seq, []).append(interval)

    total = 0
    for partner_seq, intervals in by_seq.items():
        for locus in _merge(intervals):
            if partner_seq == seq_id and _same_locus((start, end), locus):
                continue
            total += 1
    return total + 1


def annotate_repeat_families(tandems, alignments, min_loci=3):
    """Add repeat_loci and dispersed_family to each candidate."""
    if not tandems or not alignments:
        for td in tandems:
            td["repeat_loci"] = 1
            td["dispersed_family"] = False
        return tandems

    index = LocusIndex(alignments)
    for td in tandems:
        n = count_loci(index, td["seq_id"], td["copy1_start"], td["copy1_end"])
        td["repeat_loci"] = n
        td["dispersed_family"] = n >= min_loci

    n_flagged = sum(1 for td in tandems if td["dispersed_family"])
    if n_flagged:
        logger.info(
            f"  {n_flagged} candidate(s) belong to a repeat family present at "
            f"{min_loci}+ loci (likely mobile elements or rRNA operons)"
        )
    return tandems
