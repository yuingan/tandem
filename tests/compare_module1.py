"""Summarise a Module 1 run, or compare an old run with a new one.

    python tests/compare_module1.py RUN
    python tests/compare_module1.py OLD_RUN NEW_RUN [--detail]

A run is one Module 1 output folder, or a folder of them (one per genome,
matched by folder name between the two runs). Columns are read by name, so
output from v0.3.0 works as well as from later versions. --detail lists the
old calls that are now inverted pairs, multiples of a shorter repeat unit or
gone, and the new tandems the old run never reported.
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
    _, inverted = read_tsv(os.path.join(folder, "inverted_pairs.tsv"))
    _, multiples = read_tsv(os.path.join(folder, "unit_multiples.tsv"))
    return {"new_format": "architecture" in header, "calls": calls,
            "pairs": pairs, "inverted": inverted, "multiples": multiples}


def is_inverted(row):
    return row.get("orientation", "direct") == "inverted"


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
        row["inverted"] = len(run["inverted"])
        row["multiples"] = len(run["multiples"])
    else:
        row["flank"] = sum(has_flank(r) for r in calls)
        row["no_flank"] = len(calls) - row["flank"]
        # v0.3.0 kept dispersed and inverted pairs in the main table; count
        # the ones the current rules would set aside
        row["array"] = row["diverged"] = row["multiples"] = None
        row["inverted"] = sum(is_inverted(r) for r in calls)
        row["dispersed"] = sum(int(r.get("distance", 0)) > 100 and not is_inverted(r)
                               for r in calls)
    return row


def fmt(value):
    return "-" if value is None else str(value)


def print_composition(runs):
    print(f"{'genome':<28} {'tandems':>8} {'flank rep':>10} {'no flank':>9} "
          f"{'arrays':>7} {'diverged':>9} {'dispersed':>10} {'inverted':>9} {'multiples':>10}")
    for name, run in runs.items():
        c = composition(run)
        print(f"{name:<28} {c['tandems']:>8} {c['flank']:>10} "
              f"{c['no_flank']:>9} {fmt(c['array']):>7} "
              f"{fmt(c['diverged']):>9} {fmt(c['dispersed']):>10} {fmt(c['inverted']):>9} "
              f"{fmt(c['multiples']):>10}")
    if any(not r["new_format"] for r in runs.values()):
        print("  (v0.3.0 output: 'dispersed' and 'inverted' count main-table calls: direct "
              "copies more than 100 bp apart, and reversed copies)")


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


def fate(old, new):
    """Sort old calls into kept / dispersed / inverted / multiple / gone, and
    find new calls with no old match. Multiples come back as (old, new) pairs
    so the unit can be shown."""
    kept, moved, flipped, multiple, gone, matched = [], [], [], [], [], set()
    for o in old["calls"]:
        hit = next((i for i, n in enumerate(new["calls"]) if same_place(o, n)), None)
        if hit is not None:
            kept.append(o)
            matched.add(hit)
        elif any(same_place(o, p) for p in new["pairs"]):
            moved.append(o)
        elif any(same_place(o, p) for p in new["inverted"]):
            flipped.append(o)
        else:
            m = next((p for p in new["multiples"] if same_place(o, p)), None)
            if m is not None:
                multiple.append((o, m))
            else:
                gone.append(o)
    fresh = [n for i, n in enumerate(new["calls"]) if i not in matched]
    return kept, moved, flipped, multiple, gone, fresh


def print_fate(old_runs, new_runs):
    print(f"{'genome':<28} {'old calls':>10} {'still tandem':>13} "
          f"{'now dispersed':>14} {'now inverted':>13} {'now multiple':>13} {'gone':>6} "
          f"{'new tandems':>12} {'not in old':>11}")
    for name in sorted(set(old_runs) | set(new_runs)):
        if name not in old_runs or name not in new_runs:
            print(f"{name:<28} {'(only in one run)':>40}")
            continue
        kept, moved, flipped, multiple, gone, fresh = fate(old_runs[name], new_runs[name])
        print(f"{name:<28} {len(old_runs[name]['calls']):>10} {len(kept):>13} "
              f"{len(moved):>14} {len(flipped):>13} {len(multiple):>13} {len(gone):>6} "
              f"{len(new_runs[name]['calls']):>12} {len(fresh):>11}")


GAP_BINS = [(101, 500), (501, 2000), (2001, 10000), (10001, None)]


def print_gaps(runs):
    """How far apart the dispersed pairs are. The 100 bp tandem cut-off sits
    at the bottom of this; the first column is what a looser one would add."""
    labels = [f"{a}-{b}" if b else f">{a - 1}" for a, b in GAP_BINS]
    print(f"{'genome':<28} " + " ".join(f"{l:>11}" for l in labels))
    for name, run in runs.items():
        if not run["new_format"]:
            continue
        gaps = [int(p.get("distance", 0)) for p in run["pairs"]]
        counts = [sum(a <= g and (b is None or g <= b) for g in gaps) for a, b in GAP_BINS]
        print(f"{name:<28} " + " ".join(f"{c:>11}" for c in counts))


def describe(row):
    gap = row.get("distance", "?")
    extra = ""
    if "architecture" in row:
        extra = f"  {row['architecture']}  loci={row.get('repeat_loci', '?')}"
    return (f"{row['seq_id']}:{int(row['copy1_start']):,}-{int(row['copy1_end']):,} / "
            f"{int(row['copy2_start']):,}-{int(row['copy2_end']):,}  "
            f"{row.get('size', '?')}bp  {row.get('identity', '?')}%  gap {gap}{extra}")


def print_detail(old_runs, new_runs):
    for name in sorted(set(old_runs) & set(new_runs)):
        _, _, flipped, multiple, gone, fresh = fate(old_runs[name], new_runs[name])
        if not flipped and not multiple and not gone and not fresh:
            continue
        print(name)
        for row in flipped:
            print(f"  inverted  {describe(row)}")
        for row, m in multiple:
            print(f"  multiple  {describe(row)}  ({m.get('unit', '?')} bp unit)")
        for row in gone:
            print(f"  gone      {describe(row)}")
        for row in fresh:
            print(f"  new       {describe(row)}")


def main(args):
    detail = "--detail" in args
    args = [a for a in args if a != "--detail"]
    if len(args) == 1:
        runs = {k: load(v) for k, v in find_runs(args[0]).items()}
        print_composition(runs)
        print_microhomology(runs)
        print()
        print("Gaps between the copies of dispersed pairs (bp)")
        print_gaps(runs)
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
        print()
        print("Gaps between the copies of dispersed pairs, new run (bp)")
        print_gaps(new)
        if detail:
            print()
            print("Calls gone from the old run, and new tandems it never reported")
            print_detail(old, new)
    else:
        print(__doc__)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
