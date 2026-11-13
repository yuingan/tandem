#!/usr/bin/env python3
"""breseq's view of known tandem junctions, in the form the benchmark scorer reads.

    python tests/breseq_junctions.py evidence --tandem-dir module3_v0.4.0.dev13_benchmark \
        --species ecoli_K12_MG1655 --gd evidence.gd --targets targets.tsv
    python tests/breseq_junctions.py convert --breseq-dir runs/ecoli_K12_MG1655_day40 \
        --targets targets.tsv --out ecoli_K12_MG1655_day40/junction_quantification.tsv
    python tests/breseq_junctions.py check

evidence: one JC line per duplication, for --user-evidence-gd, at the junction
Module 3 measured for it on the day with the most junction reads. breseq
writes a tandem junction as the unit's first base on strand 1 and its last
base on strand -1, with 0 overlap; that is the form used here.

convert: picks, for each duplication, the JC line in breseq's output closest
to its junction and with the most reads. breseq's frequency is
J / (J + S), new-junction reads against the mean of the two reference sides;
both sides survive in a duplicated genome, so a duplication carried by a share
p of cells reads p/(1+p), and dosage is F/(1-F). With --mode denovo,
'called' is 1 when breseq accepted the junction (no reject= field).

check: runs the conversion on a JC line from breseq 0.39.0 (isolate M1, a
duplication in every cell).
"""

import argparse
import csv
import glob
import os
import sys

MATCH = 25      # bp, on top of any overlap breseq reports

# data/annotated.gd keeps every JC breseq weighed, user-supplied ones included
# (marked user_defined=1), accepted or not
GD_FILES = ("data/annotated.gd", "output/evidence/annotated.gd", "output/evidence/evidence.gd",
            "output/output.gd", "data/output.gd")

M1_LINE = ("JC\t55\t.\tNC_012660\t1848645\t1\tNC_012660\t2864810\t-1\t0\t"
           "alignment_overlap=0\tfrequency=5.098e-01\tnew_junction_read_count=26\t"
           "polymorphism_frequency=5.098e-01\tprediction=polymorphism\t"
           "side_1_read_count=27\tside_1_redundant=0\tside_2_read_count=23\tside_2_redundant=0")


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_jc(line):
    f = line.rstrip("\n").split("\t")
    if f[0] != "JC" or len(f) < 10:
        return None
    jc = {"seq1": f[3], "pos1": int(f[4]), "strand1": int(f[5]),
          "seq2": f[6], "pos2": int(f[7]), "strand2": int(f[8]),
          "overlap": int(number(f[9]) or 0)}
    for item in f[10:]:
        if "=" in item:
            key, value = item.split("=", 1)
            jc[key] = value
    return jc


def read_jcs(breseq_dir):
    jcs, seen = [], set()
    for name in GD_FILES:
        path = os.path.join(breseq_dir, name)
        if not os.path.isfile(path):
            continue
        with open(path) as f:
            for line in f:
                jc = parse_jc(line) if line.startswith("JC\t") else None
                if jc is None:
                    continue
                key = (jc["pos1"], jc["strand1"], jc["pos2"], jc["strand2"], jc["overlap"])
                if key not in seen:
                    seen.add(key)
                    jcs.append(jc)
    return jcs


def as_tandem(jc):
    """(start, end) if the JC joins a unit's end back to its start, else None."""
    sides = sorted([(jc["pos1"], jc["strand1"]), (jc["pos2"], jc["strand2"])],
                   key=lambda s: -s[1])
    (start, s1), (end, s2) = sides
    if jc["seq1"] != jc["seq2"] or (s1, s2) != (1, -1) or start >= end:
        return None
    return start, end


def best_match(jcs, seq_id, start, end):
    near = []
    for jc in jcs:
        span = as_tandem(jc)
        if span is None or jc["seq1"] != seq_id:
            continue
        tol = MATCH + jc["overlap"] + int(number(jc.get("alignment_overlap")) or 0)
        if abs(span[0] - start) <= tol and abs(span[1] - end) <= tol:
            near.append(jc)
    if not near:
        return None
    return max(near, key=lambda jc: (number(jc.get("new_junction_read_count")) or 0,
                                     -abs(as_tandem(jc)[0] - start) - abs(as_tandem(jc)[1] - end)))


def as_dosage(freq):
    if freq is None:
        return None
    return round(freq / (1 - freq), 4) if freq < 1 else float("inf")


def describe(jc):
    """dosage and counts from a JC line."""
    if jc is None:
        return {"dosage": 0.0, "called": 0, "junction_reads": 0, "reference_reads": "",
                "breseq_frequency": "", "breseq_prediction": "", "breseq_reject": "none found",
                "breseq_start": "", "breseq_end": ""}
    freq = number(jc.get("frequency", jc.get("polymorphism_frequency")))
    sides = [number(jc.get(f"side_{i}_read_count")) for i in (1, 2)
             if jc.get(f"side_{i}_redundant", "0") != "1"]
    sides = [s for s in sides if s is not None]
    start, end = as_tandem(jc)
    dosage = as_dosage(freq)
    return {
        "dosage": dosage if dosage != float("inf") else None,
        # newer breseq versions give an interval for the frequency
        "dosage_low": as_dosage(number(jc.get("frequency_lower"))),
        "dosage_high": as_dosage(number(jc.get("frequency_upper"))),
        "called": 0 if jc.get("reject") else 1,
        "junction_reads": int(number(jc.get("new_junction_read_count")) or 0),
        "reference_reads": round(sum(sides) / len(sides), 1) if sides else "",
        "breseq_frequency": "" if freq is None else freq,
        "breseq_prediction": jc.get("prediction", ""),
        "breseq_reject": jc.get("reject", ""),
        "breseq_start": start,
        "breseq_end": end,
    }


def cmd_evidence(args):
    best = {}
    for path in sorted(glob.glob(os.path.join(args.tandem_dir, f"{args.species}_day*",
                                              "junction_quantification.tsv"))):
        with open(path) as f:
            rows = csv.DictReader([line for line in f if not line.startswith("#")],
                                  delimiter="\t")
            for r in rows:
                reads = int(r["junction_reads"] or 0)
                if r["name"] not in best or reads > best[r["name"]][0]:
                    best[r["name"]] = (reads, r["seq_id"], int(r["dup_start"]), int(r["dup_end"]))
    if not best:
        sys.exit(f"No Module 3 results for {args.species} under {args.tandem_dir}")
    with open(args.gd, "w") as gd, open(args.targets, "w") as out:
        gd.write("#=GENOME_DIFF\t1.0\n")
        out.write("name\tseq_id\tstart\tend\tjunction_reads_in_tandem\n")
        for i, (name, (reads, seq_id, start, end)) in enumerate(sorted(best.items()), 1):
            gd.write(f"JC\t{i}\t.\t{seq_id}\t{start}\t1\t{seq_id}\t{end}\t-1\t0\n")
            out.write(f"{name}\t{seq_id}\t{start}\t{end}\t{reads}\n")
    print(f"{len(best)} junction(s) for {args.species} -> {args.gd}")


def cmd_convert(args):
    with open(args.targets) as f:
        targets = list(csv.DictReader(f, delimiter="\t"))
    if not any(os.path.isfile(os.path.join(args.breseq_dir, name)) for name in GD_FILES):
        sys.exit(f"No breseq output under {args.breseq_dir}")
    jcs = read_jcs(args.breseq_dir)
    if not jcs and args.mode == "given":
        print(f"WARNING: no JC lines under {args.breseq_dir}, though junctions were given",
              file=sys.stderr)
    columns = ["name", "dosage", "dosage_from", "dosage_low", "dosage_high", "called",
               "junction_reads", "reference_reads", "breseq_frequency", "breseq_prediction",
               "breseq_reject", "breseq_start", "breseq_end"]
    missing = 0
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as f:
        f.write("# breseq junction frequency F converted to dosage F/(1-F)\n")
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for t in targets:
            jc = best_match(jcs, t["seq_id"], int(t["start"]), int(t["end"]))
            missing += jc is None
            row = describe(jc)
            row.update(name=t["name"], dosage_from="breseq")
            if args.mode == "given":
                row["called"] = None    # told where to look, so not a discovery
            w.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in columns})
    note = f", {missing} with no JC line" if missing else ""
    print(f"{len(targets)} junction(s) from {len(jcs)} JC line(s){note} -> {args.out}")


def cmd_check(_):
    jc = parse_jc(M1_LINE)
    found = best_match([jc], "NC_012660", 1848645, 2864810)
    row = describe(found)
    ok = (found is not None and row["called"] == 1 and row["junction_reads"] == 26
          and row["reference_reads"] == 25.0 and abs(row["dosage"] - 1.0399) < 1e-3)
    print(f"  {'PASS' if ok else 'FAIL'}  breseq JC line for M1 converted          "
          f"F {row['breseq_frequency']} -> dosage {row['dosage']}")
    far = best_match([jc], "NC_012660", 1848645, 2864900)
    ok2 = far is None
    print(f"  {'PASS' if ok2 else 'FAIL'}  a junction 90 bp away is not matched")

    # as a junction passed with --user-evidence-gd comes back, with an interval
    given = parse_jc("JC\t7\t.\tchr\t5001\t1\tchr\t9000\t-1\t0\tfrequency=0.2\t"
                     "frequency_lower=0.15\tfrequency_upper=0.25\tnew_junction_read_count=20\t"
                     "side_1_read_count=80\tside_2_read_count=80\treject=FREQUENCY_CUTOFF\t"
                     "user_defined=1")
    row = describe(best_match([given], "chr", 5001, 9000))
    ok3 = (row["dosage"] == 0.25 and row["dosage_low"] == 0.1765
           and row["dosage_high"] == 0.3333 and row["called"] == 0)
    print(f"  {'PASS' if ok3 else 'FAIL'}  a given junction and its interval converted  "
          f"{row['dosage']} ({row['dosage_low']}-{row['dosage_high']})")
    return 0 if ok and ok2 and ok3 else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("evidence")
    e.add_argument("--tandem-dir", required=True)
    e.add_argument("--species", required=True)
    e.add_argument("--gd", required=True)
    e.add_argument("--targets", required=True)
    c = sub.add_parser("convert")
    c.add_argument("--breseq-dir", required=True)
    c.add_argument("--targets", required=True)
    c.add_argument("--out", required=True)
    c.add_argument("--mode", choices=("given", "denovo"), default="denovo",
                   help="given: breseq was passed the junctions, so 'called' is left empty")
    sub.add_parser("check")
    args = ap.parse_args()
    return {"evidence": cmd_evidence, "convert": cmd_convert, "check": cmd_check}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
