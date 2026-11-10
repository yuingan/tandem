"""Summarise a Module 1 run, or compare an old run with a new one.

    python tests/compare_module1.py RUN
    python tests/compare_module1.py OLD_RUN NEW_RUN

A run is one Module 1 output folder, or a folder of them (one per genome,
matched by folder name between the two runs). Columns are read by name, so
output from v0.3.0 works as well as from later versions.
"""

import csv
import os
import sys
from collections import Counter

SLOP = 500   # bp of coordinate wobble when matching an old call to a new one


def find_runs(path):
    """{name: folder} for one output folder or a folder of them."""
    if os.path.exists(os.path.join(path, "tandem_duplications.tsv")):
        return {os.path.basename(os.path.normpath(path)): path}
    runs = {}
    for name in sorted(os.listdir(path)):
        sub = os.path.join(path, name)
        if os.path.exists(os.path.join(sub, "tandem_duplications.tsv")):
            runs[name] = sub
    return runs


def read_tsv(path):
    if not os.path.exists(path):
        return [], []
    with open(path) as f:
        lines = [line for line in f if not line.startswith("#")]
    if not lines:
        return [], []
    reader = csv.DictReader(lines, delimiter="\t")
    rows = [r for r in reader if r.get("seq_id")]
    return reader.fieldnames or [], rows


def load(folder):
    header, calls = read_tsv(os.path.join(folder, "tandem_duplications.tsv"))
    _, pairs = read_tsv(os.path.join(folder, "dispersed_pairs.tsv"))
    return {"new_format": "architecture" in header, "calls": calls, "pairs": pairs}


def true(value):
    return str(value).strip() == "True"


def has_flank(row):
    # renamed from is_hr_signature after v0.3.0
    return true(row.get("flanking_repeat", row.get("is_hr_signature")))


def same_place(a, b):
    return (a["seq_id"] == b["seq_id"]
            and abs(int(a["copy1_start"]) - int(b["copy1_start"])) <= SLOP
            and abs(int(a["copy2_start"]) - int(b["copy2_start"])) <= SLOP)


def composition(run):
    """Counts per category; for new output flank + no flank + arrays = tandems."""
    calls = run["calls"]
    row = {"tandems": len(calls)}
    if run["new_format"]:
        arch = Counter(r.get("architecture") for r in calls)
        row["flank"] = arch["flanking_repeat"]
        row["no_flank"] = arch["no_flanking_repeat"]
        row["array"] = arch["array"]
        row["diverged"] = sum(true(r.get("copies_diverged")) for r in calls)
        row["dispersed"] = len(run["pairs"])
    else:
        row["flank"] = sum(has_flank(r) for r in calls)
        row["no_flank"] = len(calls) - row["flank"]
        # v0.3.0 kept dispersed pairs in the main table; count the ones that
        # the current gap rule would set aside
        row["array"] = row["diverged"] = None
        row["dispersed"] = sum(int(r.get("distance", 0)) > 100 for r in calls)
    return row


def fmt(value):
    return "-" if value is None else str(value)


def print_composition(runs):
    print(f"{'genome':<28} {'tandems':>8} {'flank rep':>10} {'no flank':>9} "
          f"{'arrays':>7} {'diverged':>9} {'dispersed':>10}")
    for name, run in runs.items():
        c = composition(run)
        print(f"{name:<28} {c['tandems']:>8} {c['flank']:>10} "
              f"{c['no_flank']:>9} {fmt(c['array']):>7} "
              f"{fmt(c['diverged']):>9} {fmt(c['dispersed']):>10}")
    if any(not r["new_format"] for r in runs.values()):
        print("  (v0.3.0 output: 'dispersed' counts main-table calls with a gap over 100 bp)")


def print_microhomology(runs):
    counts = Counter()
    for run in runs.values():
        for r in run["calls"]:
            if not has_flank(r) and r.get("architecture", "no_flanking_repeat") == "no_flanking_repeat":
                counts[r.get("microhomology_bp", "?")] += 1
    if counts:
        top = sorted(counts.items(), key=lambda kv: -kv[1])[:10]
        print("Microhomology at junctions without a flanking repeat: "
              + ", ".join(f"{mh} bp x{n}" for mh, n in top))


def print_fate(old_runs, new_runs):
    print(f"{'genome':<28} {'old calls':>10} {'still tandem':>13} "
          f"{'now dispersed':>14} {'gone':>6} {'new only':>9}")
    for name in sorted(set(old_runs) | set(new_runs)):
        if name not in old_runs or name not in new_runs:
            print(f"{name:<28} {'(only in one run)':>40}")
            continue
        old, new = old_runs[name], new_runs[name]
        kept = moved = gone = 0
        matched = set()
        for o in old["calls"]:
            hit = next((i for i, n in enumerate(new["calls"]) if same_place(o, n)), None)
            if hit is not None:
                kept += 1
                matched.add(hit)
            elif any(same_place(o, p) for p in new["pairs"]):
                moved += 1
            else:
                gone += 1
        new_only = len(new["calls"]) - len(matched)
        print(f"{name:<28} {len(old['calls']):>10} {kept:>13} {moved:>14} "
              f"{gone:>6} {new_only:>9}")


def main(args):
    if len(args) == 1:
        runs = {k: load(v) for k, v in find_runs(args[0]).items()}
        print_composition(runs)
        print_microhomology(runs)
    elif len(args) == 2:
        old = {k: load(v) for k, v in find_runs(args[0]).items()}
        new = {k: load(v) for k, v in find_runs(args[1]).items()}
        print("What happened to the old calls")
        print_fate(old, new)
        print()
        print("Old run")
        print_composition(old)
        print()
        print("New run")
        print_composition(new)
    else:
        print(__doc__)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
