#!/usr/bin/env python3
"""How often does the flanking-repeat search of Modules 2-3 fire at a junction without a repeat?

Random tandem duplications are drawn in the longest sequence of a reference: both ends in
sequence found once in the genome, and the 100 bp around the two ends unrelated to each
other. The function that Modules 2 and 3 use to describe a junction is then called at the
exact coordinates. No reads are needed: the description is made from the reference and
the coordinates alone. Any flanking repeat reported is therefore a repeat that lies
within the search windows but not at the junction.

Usage: random_junctions.py REFERENCE.fna N_PER_SIZE_CLASS [SEED] > out.tsv
Needs the installed tandem package (and Biopython, which it already requires).
"""
import random
import sys
import time

from Bio.Align import PairwiseAligner
from tandem.junction import classify_single_copy_junction

COMP = str.maketrans("ACGT", "TGCA")
CLASSES = [("1.5-10 kb", 1500, 10000), ("10-50 kb", 10000, 50000), ("50 kb-1 Mb", 50000, 1000000)]


def read_longest(path):
    recs, cur, seq = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    recs[cur] = "".join(seq).upper()
                cur, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if cur is not None:
        recs[cur] = "".join(seq).upper()
    name = max(recs, key=lambda k: len(recs[k]))
    return name, recs[name]


def unique_here(G, RC, pos, span=300, k=25, step=25):
    for i in range(pos - span, pos + span - k, step):
        kmer = G[i:i + k]
        if "N" in kmer or G.count(kmer) != 1 or RC.count(kmer) != 0:
            return False
    return True


def main(ref, n, seed=1):
    rng = random.Random(seed)
    name, G = read_longest(ref)
    RC = G.translate(COMP)[::-1]
    L = len(G)
    near = PairwiseAligner()
    near.mode = "local"
    near.match_score, near.mismatch_score = 2, -3
    near.open_gap_score, near.extend_gap_score = -5, -2
    print("size_class\tseq_id\tstart\tend\tsize\tflanking_repeat\trepeat_len\trepeat_identity\t"
          "geometry\tmicrohomology")
    t0 = time.time()
    for label, lo, hi in CLASSES:
        hi = min(hi, L // 3)                       # small genomes: keep the unit under a third
        if hi <= lo:
            print(f"# {label}: skipped, genome of {L:,} bp", file=sys.stderr)
            continue
        done = 0
        while done < n:
            size = int(lo * (hi / lo) ** rng.random())          # log-uniform within the class
            s = rng.randrange(20000, L - size - 20000)
            e = s + size
            if not (unique_here(G, RC, s) and unique_here(G, RC, e)):
                continue
            # the two ends must be unrelated: no local alignment worth 15 matching bases
            # between the 100 bp around the start and the 100 bp around the end
            if near.score(G[s - 50:s + 50], G[e - 50:e + 50]) >= 30:
                continue
            r = classify_single_copy_junction(G, s, e)
            print(f"{label}\t{name}\t{s + 1}\t{e}\t{size}\t{r['flanking_repeat']}\t{r['repeat_len']}\t"
                  f"{r['repeat_identity']}\t{r['repeat_geometry']}\t{r['microhomology_bp']}", flush=True)
            done += 1
        print(f"# {label}: {n} junctions, {time.time() - t0:.0f} s", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 1)
