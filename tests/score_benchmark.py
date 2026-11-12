"""Score a Module 1 run against the benchmark truth set.

    python tests/score_benchmark.py tests/data/benchmark_truth.json OUTDIR
"""

import csv
import json
import os
import sys

SLOP = 500   # bp of coordinate wobble allowed when matching a call to truth


def load_calls(outdir):
    path = os.path.join(outdir, "tandem_duplications.tsv")
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def count_dispersed(outdir):
    path = os.path.join(outdir, "dispersed_pairs.tsv")
    if not os.path.exists(path):
        return 0
    with open(path) as f:
        return sum(1 for _ in f) - 1


def matches(call, truth):
    start = int(call["copy1_start"])
    return abs(start - truth["start"]) <= SLOP


def main(truth_path, outdir):
    truth = json.load(open(truth_path))
    calls = load_calls(outdir)

    hits = []
    unmatched = list(calls)
    for t in truth:
        found = next((c for c in unmatched if matches(c, t)), None)
        if found:
            unmatched.remove(found)
        hits.append((t, found))

    print(f"{outdir}")
    print(f"  {sum(1 for _, c in hits if c)}/{len(truth)} duplications recovered")
    for t, c in hits:
        if c:
            size_ok = abs(int(c["size"]) - t["size"]) <= SLOP
            copies_ok = round(float(c.get("copy_number", 2))) == t["copies"]
            flags = "" if (size_ok and copies_ok) else "  <- size/copies wrong"
            print(f"    {t['size']:>6,}bp x{t['copies']} {t['identity']:>5}%  "
                  f"called {c['size']:>6}bp x{c.get('copy_number','?')} "
                  f"{c.get('architecture','?'):<20} div={c.get('copies_diverged','?'):<5} loci={c.get('repeat_loci','?')}{flags}")
        else:
            print(f"    {t['size']:>6,}bp x{t['copies']} {t['identity']:>5}%  MISSED")

    print(f"  {len(unmatched)} call(s) not in the truth set")
    for c in unmatched[:10]:
        print(f"    {c['seq_id']}:{c['copy1_start']} {c['size']}bp "
              f"loci={c.get('repeat_loci','?')} family={c.get('dispersed_family','?')}")
    print(f"  {count_dispersed(outdir)} dispersed pair(s) set aside")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
