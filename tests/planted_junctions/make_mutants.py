#!/usr/bin/env python3
"""Plant fresh tandem duplications with known junctions in a genome.

Ten mutant genomes are written, each a copy of the reference with one tandem
duplication of 1.5-48 kb in its longest sequence. The junctions follow the design of
the artificial-population benchmark:

  m01      no microhomology, two bases inserted at the junction (chosen so that they
           cannot be read as part of either copy)
  m02      no microhomology
  m03-m05  1, 2 and 4 bp of microhomology already present in the genome
  m06-m08  7, 12 and 18 bp of microhomology (the bases before the duplicated unit are
           replaced by the unit's last bases)
  m09-m10  a duplication formed by a crossover between direct repeats of 60 and 80 bp
           at 95% identity (the second repeat is planted after the unit)

Both ends of every duplicated unit lie in sequence found once in the genome.

Usage: make_mutants.py REFERENCE.fna OUT_DIR [SEED]
Writes OUT_DIR/m01.fna ... m10.fna and OUT_DIR/truth.tsv (1-based, inclusive
coordinates in the mutant genome). Standard library only.
"""
import os
import random
import sys

COMP = str.maketrans("ACGT", "TGCA")

DESIGNS = [  # id, kind, microhomology or repeat length, approximate unit size
    ("m01", "insertion", 0, 48000),
    ("m02", "mh", 0, 30000),
    ("m03", "mh", 1, 12000),
    ("m04", "mh", 2, 5000),
    ("m05", "mh", 4, 1500),
    ("m06", "mh", 7, 20000),
    ("m07", "mh", 12, 8000),
    ("m08", "mh", 18, 3000),
    ("m09", "repeat", 60, 25000),
    ("m10", "repeat", 80, 10000),
]


def read_fasta(path):
    records, name, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq).upper()))
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        records.append((name, "".join(seq).upper()))
    return records


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


def unique_here(G, RC, pos, span=300, k=25, step=25):
    """True if the 25-mers sampled within span of pos occur once in the sequence."""
    for i in range(pos - span, pos + span - k, step):
        kmer = G[i:i + k]
        if "N" in kmer or G.count(kmer) != 1 or RC.count(kmer) != 0:
            return False
    return True


def mutate(seq, n, rng):
    """Copy of seq with n substitutions spread evenly along it."""
    s = list(seq)
    for j in range(n):
        i = int((j + 0.5) * len(s) / n)
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def main(ref, out, seed=1):
    rng = random.Random(seed)
    records = read_fasta(ref)
    main_i = max(range(len(records)), key=lambda i: len(records[i][1]))
    name, G = records[main_i]
    RC = G.translate(COMP)[::-1]
    L = len(G)
    margin = min(200_000, L // 10)
    os.makedirs(out, exist_ok=True)
    rows = []
    for mid, kind, k, size in DESIGNS:
        for _ in range(100000):
            s = rng.randrange(margin, L - margin - size)
            e = s + size + rng.randrange(0, 400)
            if not (unique_here(G, RC, s) and unique_here(G, RC, e)):
                continue
            if kind == "repeat":
                r = k
                mismatches = round(r * 0.05)
                R1 = G[s - r:s]
                R2 = mutate(R1, mismatches, rng)
                Rh = R2[:r // 2] + R1[r // 2:]          # crossover in the middle of the repeat
                D = G[s:e]
                M = G[:s - r] + R1 + D + Rh + D + R2 + G[e + r:]
                rows.append(dict(id=mid, seq_id=name, design=kind, microhomology="", inserted=0,
                                 repeat_len=r, repeat_identity=round(100 * (r - mismatches) / r, 1),
                                 copy1_start=s - r + 1, copy1_end=e, copy2_start=e + 1,
                                 copy2_end=2 * e - s + r, unit=e - s + r))
                break
            if kind == "mh" and 1 <= k <= 4:
                # microhomology of exactly k bp as the genome has it: move the end until the
                # k bases before it equal the k bases before the start, and no further
                found = None
                for e2 in range(e, e + 4000):
                    if (G[e2 - k:e2] == G[s - k:s] and G[e2 - k - 1] != G[s - k - 1]
                            and G[e2] != G[s]):
                        found = e2
                        break
                if found is None or not unique_here(G, RC, found):
                    continue
                e = found
                left = G[:s]
            elif k == 0:
                if G[s - 1] == G[e - 1] or G[s] == G[e]:
                    continue
                left = G[:s]
            else:
                if G[s - k - 1] == G[e - k - 1] or G[s] == G[e]:
                    continue
                left = G[:s - k] + G[e - k:e]
            D = G[s:e]
            ins = ""
            if kind == "insertion":
                # two bases that cannot be read as part of either copy: the first differs
                # from the base after the unit, the second from the base before it
                ins = (rng.choice([b for b in "ACGT" if b != G[e]])
                       + rng.choice([b for b in "ACGT" if b != G[s - 1]]))
            M = left + D + ins + D + G[e:]
            rows.append(dict(id=mid, seq_id=name, design=kind, microhomology=k, inserted=len(ins),
                             repeat_len=0, repeat_identity="", copy1_start=s + 1, copy1_end=e,
                             copy2_start=e + 1 + len(ins), copy2_end=2 * e - s + len(ins),
                             unit=e - s))
            break
        else:
            sys.exit(f"{mid}: no suitable position found in {name}")
        mutant = list(records)
        mutant[main_i] = (name, M)
        write_fasta(os.path.join(out, f"{mid}.fna"), mutant)
    keys = list(rows[0].keys())
    with open(os.path.join(out, "truth.tsv"), "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in keys) + "\n")
    print(f"{ref}: {len(rows)} mutant genomes in {out}/ ({name}, {L:,} bp)")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 1)
