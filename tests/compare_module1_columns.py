"""Compare two Module 1 runs call by call and column by column.

    python tests/compare_module1_columns.py OLD_RUN NEW_RUN

A run is one Module 1 output folder, or a folder of them (one per genome,
matched by folder name). Calls are matched by their coordinates in each of the
five tables. For every call present in both runs the columns that differ are
listed, and calls present in one run only are listed too.

Written for 0.4.0.dev20, which changed the flanking-repeat search of pairs
with a unit under 400 bp and nothing else. The last line says whether the
differences are of that kind only: the same calls in both runs, and changes
confined to the flanking-repeat columns of units under 400 bp.
"""

import csv
import os
import sys

TABLES = ("tandem_duplications.tsv", "dispersed_pairs.tsv", "inverted_pairs.tsv",
          "unit_multiples.tsv", "short_unit_arrays.tsv")
KEY = ("seq_id", "copy1_start", "copy1_end", "copy2_start", "copy2_end")
REPEAT_COLUMNS = {"flanking_repeat", "repeat_len", "repeat_identity", "repeat_geometry",
                  "architecture", "classification_note"}


def find_runs(path):
    """{name: folder} for one output folder or a folder of them."""
    if os.path.exists(os.path.join(path, TABLES[0])):
        return {os.path.basename(os.path.normpath(path)): path}
    runs = {}
    for name in sorted(os.listdir(path)):
        sub = os.path.join(path, name)
        if os.path.exists(os.path.join(sub, TABLES[0])):
            runs[name] = sub
    return runs


def read_table(path):
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        lines = [line for line in f if not line.startswith("#")]
    if not lines:
        return {}
    rows = {}
    for r in csv.DictReader(lines, delimiter="\t"):
        if r.get("seq_id"):
            rows[tuple(r.get(k, "") for k in KEY)] = r
    return rows


def is_array(row):
    if row.get("architecture"):
        return row["architecture"] == "array"
    return float(row.get("copy_number") or 2) > 2


def table1(calls):
    """All, pairs, arrays and calls with copies at 99% identity or above."""
    calls = list(calls)
    arrays = sum(is_array(r) for r in calls)
    recent = sum(float(r.get("identity") or 0) >= 99 for r in calls)
    return len(calls), len(calls) - arrays, arrays, recent


def flanked_pairs(calls):
    return sum(r.get("flanking_repeat") == "True" and not is_array(r) for r in calls)


def main(old_path, new_path):
    old_runs, new_runs = find_runs(old_path), find_runs(new_path)
    names = [n for n in new_runs if n in old_runs]
    if not names:
        sys.exit(f"no genome folder with {TABLES[0]} in both {old_path} and {new_path}")
    for n in sorted(set(old_runs) ^ set(new_runs)):
        print(f"{n}: in one run only, not compared")

    unexpected = 0
    changed_calls = 0
    print(f"{'genome':<24}{'all':>5}{'pairs':>7}{'arrays':>8}{'>=99%':>7}"
          f"{'pairs with a flanking repeat, old -> new':>44}")
    details = []
    for name in names:
        new_calls = read_table(os.path.join(new_runs[name], TABLES[0]))
        old_calls = read_table(os.path.join(old_runs[name], TABLES[0]))
        a, p, r, c = table1(new_calls.values())
        print(f"{name:<24}{a:>5}{p:>7}{r:>8}{c:>7}"
              f"{flanked_pairs(old_calls.values()):>38} -> {flanked_pairs(new_calls.values())}")
        if table1(old_calls.values()) != (a, p, r, c):
            details.append(f"{name}: counts were {table1(old_calls.values())} in the old run")
            unexpected += 1
        for table in TABLES:
            old = read_table(os.path.join(old_runs[name], table))
            new = read_table(os.path.join(new_runs[name], table))
            for key in sorted(set(old) - set(new)):
                details.append(f"{name} {table}: {key[1]}-{key[4]} in the old run only")
                unexpected += 1
            for key in sorted(set(new) - set(old)):
                details.append(f"{name} {table}: {key[1]}-{key[4]} in the new run only")
                unexpected += 1
            for key in sorted(set(old) & set(new), key=lambda k: int(k[1])):
                o, w = old[key], new[key]
                cols = [k for k in w if k in o and o[k] != w[k]]
                cols += [k for k in set(o) ^ set(w)]
                if not cols:
                    continue
                changed_calls += 1
                size = int(float(w.get("size") or 0))
                fine = set(cols) <= REPEAT_COLUMNS and size < 400
                unexpected += not fine
                what = "; ".join(f"{k} {o.get(k, '-')} -> {w.get(k, '-')}" for k in cols
                                 if k in ("flanking_repeat", "repeat_len", "architecture")
                                 or k not in REPEAT_COLUMNS)
                details.append(f"{name} {table.split('.')[0]}: {key[1]}, unit {size} bp, "
                               f"x{w.get('copy_number', '?')}: {what}"
                               + ("" if fine else "   <-- not expected"))
    print()
    for line in details:
        print("  " + line)
    print()
    if unexpected:
        print(f"UNEXPECTED: {unexpected} difference(s) are not of the expected kind (marked above).")
        return 1
    print(f"As expected: the same calls in both runs; {changed_calls} call(s) differ, all with a unit "
          "under 400 bp and only in the flanking-repeat columns.")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1], sys.argv[2]))
