"""Simulate Illumina-like reads from a genome, with a duplication at a set
frequency, so junction counts can be checked against a number we chose.

    python tests/simulate_reads.py REF.fna OUT_PREFIX --start S --end E \
        --frequency 0.3 --copies 2 --depth 100

Coordinates are 1-based inclusive over the reference, and describe the unit
that gets duplicated in the fraction of the population given by --frequency.
"""

import argparse
import gzip
import math
import random
import sys



COMPLEMENT = str.maketrans("ACGT", "TGCA")


def load_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def rc(seq):
    return seq.translate(COMPLEMENT)[::-1]


def sequence_error(seq, rate):
    """Place errors by drawing how many, not by rolling a die per base."""
    if rate <= 0:
        return seq
    n = _errors_per_read(rate, len(seq))
    if n == 0:
        return seq
    bases = list(seq)
    for i in random.sample(range(len(bases)), n):
        bases[i] = random.choice([b for b in "ACGT" if b != bases[i]])
    return "".join(bases)


def _errors_per_read(rate, length):
    """Poisson draw, which is close enough to binomial at these rates."""
    lam = rate * length
    n, p, target = 0, math.exp(-lam), random.random()
    cumulative = p
    while cumulative < target and n < length:
        n += 1
        p *= lam / n
        cumulative += p
    return n


def emit(fh, name, seq):
    fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reference")
    ap.add_argument("prefix")
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--frequency", type=float, default=1.0)
    ap.add_argument("--copies", type=int, default=2)
    ap.add_argument("--depth", type=float, default=100)
    ap.add_argument("--read-length", type=int, default=150)
    ap.add_argument("--insert", type=int, default=350)
    ap.add_argument("--error", type=float, default=0.001)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--gzip", action="store_true")
    args = ap.parse_args()
    random.seed(args.seed)

    ref = load_fasta(args.reference)
    unit = ref[args.start - 1:args.end]

    wild = ref
    mutant = ref[:args.end] + unit * (args.copies - 1) + ref[args.end:]

    n_wild = int(len(wild) * args.depth * (1 - args.frequency) / args.read_length / 2)
    n_mut = int(len(mutant) * args.depth * args.frequency / args.read_length / 2)

    opener = gzip.open if args.gzip else open
    suffix = ".fq.gz" if args.gzip else ".fq"
    r1 = opener(f"{args.prefix}_R1{suffix}", "wt")
    r2 = opener(f"{args.prefix}_R2{suffix}", "wt")

    idx = 0
    for template, n in ((wild, n_wild), (mutant, n_mut)):
        span = len(template) - args.insert
        if span <= 0 or n <= 0:
            continue
        for _ in range(n):
            pos = random.randrange(span)
            frag = template[pos:pos + args.insert]
            fwd = sequence_error(frag[:args.read_length], args.error)
            rev = sequence_error(rc(frag[-args.read_length:]), args.error)
            emit(r1, f"r{idx}/1", fwd)
            emit(r2, f"r{idx}/2", rev)
            idx += 1

    r1.close()
    r2.close()
    print(f"{args.prefix}_R1{suffix}  {idx:,} pairs  "
          f"unit {args.end - args.start + 1:,}bp  f={args.frequency}  n={args.copies}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
