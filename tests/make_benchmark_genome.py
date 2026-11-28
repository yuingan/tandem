"""Build one genome with a known answer and a realistic repeat background.

Toy genomes with a single repeat in them do not tell you much: they are the
easiest possible case for detection and the hardest for NUCmer's simplify
step. This builds something closer to a bacterial chromosome, with rRNA
operons, two IS families and duplications spanning three orders of magnitude,
then writes the truth set next to it.

    python tests/make_benchmark_genome.py tests/data
"""

import json
import os
import random
import sys

random.seed(2718)

BACKBONE = 1_500_000


def rand_seq(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def mutate(seq, identity):
    """Return a copy of seq at roughly the given percent identity."""
    if identity >= 100:
        return seq
    bases = list(seq)
    n_changes = int(len(bases) * (100 - identity) / 100)
    for pos in random.sample(range(len(bases)), n_changes):
        bases[pos] = random.choice([b for b in "ACGT" if b != bases[pos]])
    return "".join(bases)


def build():
    """Assemble the chromosome from a list of inserts placed on a backbone."""
    backbone = rand_seq(BACKBONE)

    rrna = rand_seq(5000)
    is_a = rand_seq(1300)
    is_b = rand_seq(900)

    # (sort key, sequence, truth record or None)
    inserts = []

    # rRNA operons: seven near-identical copies, scattered. Real repeats,
    # not duplications.
    for i in range(7):
        pos = 60_000 + i * 190_000
        inserts.append((pos, mutate(rrna, 99.5), None))

    # Two IS families at different copy numbers. Real chromosomes carry
    # plenty of same-family copies a few kilobases apart, and those are the
    # pairs that get mistaken for duplications, so some sit well inside the
    # default search radius.
    for i in range(15):
        inserts.append((25_000 + i * 95_000, is_a, None))
    inserts.append((33_000, is_a, None))       # 8 kb from the first copy
    inserts.append((218_000, is_a, None))      # 3 kb from another
    for i in range(6):
        inserts.append((140_000 + i * 210_000, mutate(is_b, 99.0), None))
    inserts.append((160_000, mutate(is_b, 99.0), None))   # 20 kb apart

    # The duplications we expect to get back. Size, identity, copies.
    wanted = [
        (110_000, 800, 100.0, 2),
        (330_000, 2_000, 100.0, 2),
        (555_000, 5_000, 100.0, 2),
        (780_000, 12_000, 99.8, 2),
        (1_010_000, 30_000, 100.0, 2),
        (1_260_000, 4_000, 100.0, 3),
        (1_400_000, 3_000, 97.0, 2),
        (1_450_000, 6_000, 92.0, 2),
    ]
    for pos, size, identity, copies in wanted:
        unit = rand_seq(size)
        block = unit + "".join(mutate(unit, identity) for _ in range(copies - 1))
        inserts.append((pos, block, {
            "size": size, "identity": identity, "copies": copies,
        }))

    # Splice everything in, tracking where each insert landed.
    inserts.sort(key=lambda x: x[0])
    seq = []
    truth = []
    cursor = 0
    offset = 0
    for pos, block, record in inserts:
        seq.append(backbone[cursor:pos])
        start = pos + offset + 1          # 1-based
        seq.append(block)
        if record:
            record["start"] = start
            record["end"] = start + len(block) - 1
            truth.append(record)
        offset += len(block)
        cursor = pos
    seq.append(backbone[cursor:])

    return "".join(seq), truth


def main(outdir):
    os.makedirs(outdir, exist_ok=True)
    seq, truth = build()

    fasta = os.path.join(outdir, "benchmark.fna")
    with open(fasta, "w") as f:
        f.write(">benchmark\n")
        for i in range(0, len(seq), 70):
            f.write(seq[i:i + 70] + "\n")

    with open(os.path.join(outdir, "benchmark_truth.json"), "w") as f:
        json.dump(truth, f, indent=2)

    print(f"{fasta}\t{len(seq):,} bp\t{len(truth)} duplications")
    for t in truth:
        print(f"  {t['start']:>9,}  {t['size']:>6,}bp  x{t['copies']}  {t['identity']}%")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "tests/data")
