"""Flanking repeat of a tandem pair in a genome sequence, on planted sequences.

    python tests/check_flanking_repeat.py [WORKDIR DATADIR]

The search compares windows at the outer ends of the two copies. One of them
reaches back into the body of a copy, and up to 0.4.0.dev19 it reached 200 bp
whatever the unit length. For a unit under 400 bp that is more than half of
the unit, so the two windows held the same stretch of the two copies, the
copies were aligned with each other, and a pair with a unit under about
330 bp was given a flanking repeat of 200 - unit/2 bp whether or not it had
one. These checks plant pairs with and without a repeat and read the answer.

With WORKDIR and DATADIR (as run_tests.sh gives them) the suite's 150 bp pair
is also run through the command with --min-size 100.
"""

import csv
import os
import random
import subprocess
import sys

from tandem.junction import classify_junction


def rnd(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def pair(unit, before=None, after=None, second=None):
    """before + unit + (second or unit) + after, and the 0-based copy coordinates."""
    before = rnd(3000) if before is None else before
    after = rnd(3000) if after is None else after
    second = unit if second is None else second
    a = len(before)
    seq = before + unit + second + after
    return seq, a, a + len(unit), a + len(unit), a + len(unit) + len(second)


def main(work=None, data=None):
    random.seed(11)
    results = []

    # Random units with random sequence on both sides: no repeat to find.
    lengths = (150, 200, 250, 320, 390, 600)
    seen = {}
    for n in lengths:
        seq, c1s, c1e, c2s, c2e = pair(rnd(n))
        r = classify_junction(seq, c1s, c1e, c2s, c2e)
        seen[n] = r["repeat_len"] if r["flanking_repeat"] else 0
    ok = not any(seen.values())
    detail = ("none at " + ", ".join(str(n) for n in lengths) + " bp" if ok else
              "; ".join(f"{n} bp unit: {v} bp" for n, v in seen.items() if v))
    results.append(("short pair without a repeat", ok, detail))

    # R-D-R-D-R with a 60 bp repeat and a 300 bp unit, in both bracketings:
    # the repeat at the start of each copy and once more after the second,
    # or at the end of each copy and once more before the first.
    rep, body = rnd(60), rnd(240)
    seq, c1s, c1e, c2s, c2e = pair(rep + body, after=rep + rnd(3000))
    first = classify_junction(seq, c1s, c1e, c2s, c2e)
    seq, c1s, c1e, c2s, c2e = pair(body + rep, before=rnd(3000) + rep)
    second = classify_junction(seq, c1s, c1e, c2s, c2e)
    ok = all(r["flanking_repeat"] and 55 <= r["repeat_len"] <= 70 for r in (first, second))
    results.append(("short pair with a 60 bp flanking repeat", ok,
                    f"{first['repeat_len']} and {second['repeat_len']} bp"))

    # Second copy 60 bp shorter than the first: the windows are set from the
    # shorter copy, so that they do not overlap on it.
    unit = rnd(300)
    seq, c1s, c1e, c2s, c2e = pair(unit, second=unit[:120] + unit[180:])
    r = classify_junction(seq, c1s, c1e, c2s, c2e)
    results.append(("copies of unequal length", not r["flanking_repeat"],
                    f"{r['repeat_len']} bp" if r["flanking_repeat"] else "none"))

    # Units of 400 bp or more were never affected: same answers as before.
    seq, c1s, c1e, c2s, c2e = pair(rnd(5000))
    plain = classify_junction(seq, c1s, c1e, c2s, c2e)
    rep, body = rnd(80), rnd(4920)
    seq, c1s, c1e, c2s, c2e = pair(rep + body, after=rep + rnd(3000))
    flanked = classify_junction(seq, c1s, c1e, c2s, c2e)
    ok = not plain["flanking_repeat"] and flanked["flanking_repeat"] and 75 <= flanked["repeat_len"] <= 90
    results.append(("5 kb pair, without and with an 80 bp repeat", ok,
                    f"none; {flanked['repeat_len']} bp"))

    if work and data:
        out = os.path.join(work, "tandem_short_150_min100")
        subprocess.run(["tandem", "-r", os.path.join(data, "tandem_short_150.fna"), "-o", out,
                        "-t", "1", "--min-size", "100"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
        path = os.path.join(out, "tandem_duplications.tsv")
        rows = []
        if os.path.exists(path):
            with open(path) as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
        ok = (len(rows) == 1 and 148 <= int(rows[0]["size"]) <= 152
              and rows[0]["flanking_repeat"] == "False")
        detail = (f"{rows[0]['size']} bp, flanking repeat {rows[0]['flanking_repeat']}"
                  if len(rows) == 1 else f"{len(rows)} call(s)")
        results.append(("150 bp pair at --min-size 100", ok, detail))

    failed = 0
    for name, ok, detail in results:
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<42} {detail}")
        failed += not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:3]))
