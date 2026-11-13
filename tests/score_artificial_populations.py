#!/usr/bin/env python3
"""Score Module 3 on the artificial populations against their ground truth.

    python tests/score_artificial_populations.py --gt-dir genome/artificial_pop_seqs \
        --results-dir module3_v0.4.0.dev12_benchmark

Reads <gt-dir>/<species>/ground_truth.json and every
<results-dir>/<species>_day<d>/junction_quantification.tsv. Output from
v0.3.0, which has a JRR column instead of dosage, is scored the same way, so
an old and a new run can be set side by side. All points go to
<results-dir>/scatter_data.tsv.

The truth is each mutant's share of the genomes in the simulated reads. The
generator gives the wild type coverage x (1 - sum of mutant frequencies),
floored at 0, and each mutant coverage x its frequency, so where the
frequencies of a day add up to more than 1 there is no wild type and each
mutant's share is its frequency over that sum. --labels-as-written scores
against the frequencies as written instead.

Per species, per mechanism and overall:

  recall      estimates above 1% among mutants truly above 1%
  MAE         mean absolute error over those same points
  slope       least-squares slope through the origin, estimate on truth,
              over every point where the mutant is present
  R2          of a straight line with an intercept, same points
  absent      the largest estimate where the mutant is not in the population
"""

import argparse
import csv
import json
import os
import sys


def number(value):
    return float(value) if value not in ("", "None", None) else None


def read_results(path):
    with open(path) as f:
        lines = [line for line in f if line.strip() and not line.startswith("#")]
    rows = list(csv.DictReader(lines, delimiter="\t"))
    column = "dosage" if rows and "dosage" in rows[0] else "JRR"
    out = {}
    for r in rows:
        estimate = number(r.get(column))
        source = r.get("dosage_from", "junction" if column == "dosage" else "JRR")
        out[r["name"]] = {
            "estimate": estimate,
            "source": source if estimate is not None else "",
            "junction": estimate if source in ("junction", "JRR") else None,
            "depth": number(r.get("depth_dosage")),
            "junction_reads": r.get("junction_reads", r.get("spanning_reads", "")),
            "reference_reads": r.get("reference_reads", ""),
        }
    return out, column


def load_points(gt_dir, results_dir, labels_as_written=False):
    points, columns, missing = [], set(), []
    for species in sorted(os.listdir(gt_dir)):
        gt_path = os.path.join(gt_dir, species, "ground_truth.json")
        if not os.path.isfile(gt_path):
            continue
        with open(gt_path) as f:
            truth = json.load(f)
        days = sorted({int(k.split("_")[1]) for t in truth for k in t["frequencies"]})
        for day in days:
            label = {t["mutant_id"]: float(t["frequencies"].get(f"day_{day}", 0)) for t in truth}
            total = sum(label.values())
            scale = 1.0 if labels_as_written else max(1.0, total)
            tsv = os.path.join(results_dir, f"{species}_day{day}", "junction_quantification.tsv")
            if not os.path.isfile(tsv):
                missing.append(f"{species}_day{day}")
                continue
            results, column = read_results(tsv)
            columns.add(column)
            for t in truth:
                r = results.get(t["mutant_id"], {})
                points.append({
                    "species": species,
                    "mutant_id": t["mutant_id"],
                    "mechanism": t.get("mechanism", ""),
                    "dup_size": t.get("dup_size", ""),
                    "day": day,
                    "label": label[t["mutant_id"]],
                    "true": round(label[t["mutant_id"]] / scale, 4),
                    "estimate": r.get("estimate"),
                    "source": r.get("source", ""),
                    "junction": r.get("junction"),
                    "depth": r.get("depth"),
                    "junction_reads": r.get("junction_reads", ""),
                    "reference_reads": r.get("reference_reads", ""),
                })
    return points, columns, missing


def summarise(points, key="estimate"):
    present = [p for p in points if p["true"] > 0.01]
    found = sum(1 for p in present if (p[key] or 0) > 0.01)
    mae = (sum(abs((p[key] or 0) - p["true"]) for p in present) / len(present)
           if present else None)

    fit = [p for p in points if p["true"] > 0 and p[key] is not None]
    xs = [p["true"] for p in fit]
    ys = [p[key] for p in fit]
    slope = r2 = None
    if len(fit) >= 2:
        sxx = sum(x * x for x in xs)
        slope = sum(x * y for x, y in zip(xs, ys)) / sxx if sxx else None
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
        vx = sum((x - mx) ** 2 for x in xs)
        vy = sum((y - my) ** 2 for y in ys)
        cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
        r2 = cov * cov / (vx * vy) if vx and vy else None

    absent = [p[key] for p in points if p["true"] == 0 and p[key] is not None]
    unmeasured = sum(1 for p in points if p[key] is None)
    return {
        "points": len(points),
        "recall": f"{found}/{len(present)}",
        "mae": mae,
        "slope": slope,
        "r2": r2,
        "absent_max": max(absent) if absent else None,
        "unmeasured": unmeasured,
    }


def fmt(value, spec):
    return "-" if value is None else format(value, spec)


def print_table(groups, key="estimate"):
    print(f"{'':<26} {'points':>7} {'recall':>9} {'MAE':>8} {'slope':>7} "
          f"{'R2':>6} {'absent max':>11} {'unmeasured':>11}")
    for name, pts in groups:
        s = summarise(pts, key)
        print(f"{name:<26} {s['points']:>7} {s['recall']:>9} {fmt(s['mae'], '.4f'):>8} "
              f"{fmt(s['slope'], '.3f'):>7} {fmt(s['r2'], '.3f'):>6} "
              f"{fmt(s['absent_max'], '.3f'):>11} {s['unmeasured']:>11}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--gt-dir", default="genome/artificial_pop_seqs")
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--labels-as-written", action="store_true",
                    help="score against the frequencies as written, not the shares simulated")
    args = ap.parse_args()

    points, columns, missing = load_points(args.gt_dir, args.results_dir,
                                           args.labels_as_written)
    if not points:
        sys.exit(f"No results found under {args.results_dir}")
    print(f"Scored column: {', '.join(sorted(columns))}")
    if args.labels_as_written:
        print("Truth: the frequencies as written in ground_truth.json")
    else:
        rescaled = sorted({(p["species"], p["day"]) for p in points if p["true"] != p["label"]})
        days = sorted({d for _, d in rescaled})
        print("Truth: each mutant's share of the simulated genomes"
              + (f" (frequencies rescaled on day {', '.join(map(str, days))}, "
                 f"where they add up to more than 1)" if days else ""))
    if missing:
        print(f"No results for: {', '.join(missing)}")
    print()

    species = sorted({p["species"] for p in points})
    mechanisms = sorted({p["mechanism"] for p in points if p["mechanism"]})
    groups = [(s, [p for p in points if p["species"] == s]) for s in species]
    groups += [(f"mechanism {m}", [p for p in points if p["mechanism"] == m]) for m in mechanisms]
    groups.append(("all", points))
    print_table(groups)

    sources = {}
    for p in points:
        if p["source"]:
            sources[p["source"]] = sources.get(p["source"], 0) + 1
    if len(sources) > 1 or any(p["depth"] is not None for p in points):
        print()
        print("All points, by estimate:")
        rows = [("as reported", "estimate"), ("junction reads only", "junction"),
                ("depth only", "depth")]
        print_table([(name, [dict(p, estimate=p[key]) for p in points]) for name, key in rows])
        print("  reported estimates: " + ", ".join(f"{n} from {s}" for s, n in sorted(sources.items())))

    print()
    print("By day (as reported):")
    print_table([(f"day {d}", [p for p in points if p["day"] == d])
                 for d in sorted({p["day"] for p in points})])

    print()
    print("recall and MAE over points truly above 1%; slope (through the origin) and R2")
    print("over every point where the mutant is present; 'absent max' is the largest")
    print("estimate where it is not; 'unmeasured' had no estimate at all.")

    out = os.path.join(args.results_dir, "scatter_data.tsv")
    fields = ["species", "mutant_id", "mechanism", "dup_size", "day", "label", "true",
              "estimate", "source", "junction_reads", "reference_reads", "depth"]
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for p in points:
            w.writerow({k: ("" if p[k] is None else p[k]) for k in fields})
    print(f"\nAll points: {out}")


if __name__ == "__main__":
    main()
