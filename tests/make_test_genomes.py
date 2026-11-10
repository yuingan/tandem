"""Build small genomes with known answers.

Each one targets a case a user will actually hit, including the three that
came back in review: an IS pair masquerading as a duplication, a three-copy
array, and a repeat below the size floor.

    python tests/make_test_genomes.py tests/data
"""

import os
import random
import sys

random.seed(11)


def rand_seq(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def write(path, name, seq):
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 70):
            f.write(seq[i:i + 70] + "\n")


def build():
    backbone = rand_seq(60000)
    unit = rand_seq(5000)
    is_element = rand_seq(1200)

    cases = {}

    # Straightforward two-copy tandem. Should be called.
    cases["tandem_2x"] = (
        backbone[:20000] + unit + unit + backbone[20000:],
        "one tandem duplication, 5000bp, gap 0",
    )

    # Three copies in a row. Reviewer saw nothing reported for this.
    cases["tandem_3x"] = (
        backbone[:20000] + unit * 3 + backbone[20000:],
        "one tandem array, 5000bp unit, 3 copies",
    )

    # Two IS copies 30kb apart in the same orientation. Not a duplication:
    # the 30kb between them is present once.
    cases["is_pair_30kb"] = (
        backbone[:10000] + is_element + backbone[10000:40000]
        + is_element + backbone[40000:],
        "dispersed pair, not a duplication",
    )

    # The same IS element at six sites. Every pairwise combination within
    # the search radius used to surface as a separate candidate.
    seq = backbone[:5000]
    for i in range(6):
        seq += is_element + backbone[5000 + i * 6000:5000 + (i + 1) * 6000]
    cases["is_family_6x"] = (seq, "repeat family at 6 loci, no duplications")

    # Tandem duplication of a diverged unit. Junction signature is gone;
    # the event is real but the mechanism is not callable.
    old = list(unit)
    for pos in random.sample(range(len(old)), int(len(old) * 0.05)):
        old[pos] = random.choice("ACGT")
    cases["tandem_diverged"] = (
        backbone[:20000] + unit + "".join(old) + backbone[20000:],
        "one tandem duplication at 95% identity, mechanism indeterminate",
    )

    # Short repeat, below the default size floor.
    short = rand_seq(150)
    cases["tandem_short_150"] = (
        backbone[:20000] + short + short + backbone[20000:],
        "150bp tandem, below --min-size 200",
    )

    return cases


def main(outdir):
    os.makedirs(outdir, exist_ok=True)
    for name, (seq, expected) in build().items():
        path = os.path.join(outdir, f"{name}.fna")
        write(path, name, seq)
        print(f"{path}\t{len(seq):,}bp\t{expected}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "tests/data")
