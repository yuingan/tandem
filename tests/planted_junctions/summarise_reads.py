#!/usr/bin/env python3
"""Totals over every reads_planted.tsv and random_junctions.tsv below a folder.
Usage: summarise_reads.py FOLDER"""
import csv
import os
import statistics
import sys

PLANTED_DOSAGE = 0.3    # the share of cells with the duplication in reads_planted.py


def main(top):
    rows, rand = [], []
    for d, _, files in sorted(os.walk(top)):
        if "reads_planted.tsv" in files:
            rows += list(csv.DictReader(open(os.path.join(d, "reads_planted.tsv")), delimiter="\t"))
        if "random_junctions.tsv" in files:
            g = os.path.basename(d)
            for r in csv.DictReader(open(os.path.join(d, "random_junctions.tsv")), delimiter="\t"):
                r["genome"] = g
                rand.append(r)
    if not rows:
        sys.exit(f"no reads_planted.tsv below {top}")
    genomes = sorted({r["genome"] for r in rows})
    scan = {}
    for r in rows:
        if r["source"] == "scan":
            scan.setdefault((r["genome"], r["id"]), []).append(r)

    def scan_hit(r):
        """The scan's call at the junction Module 3 reports, or None."""
        hit = [c for c in scan.get((r["genome"], r["id"]), [])
               if c["start"] and c["start"] == r["start"] and c["end"] == r["end"]]
        return hit[0] if hit else None

    def line(label, rs):
        m3 = [r for r in rs if r["source"] == "module3"]
        mh = [r for r in m3 if r["design"] == "microhomology"]
        rp = [r for r in m3 if r["design"] == "repeat"]
        same = other = scan_mh = 0
        for r in m3:
            calls = [c for c in scan.get((r["genome"], r["id"]), []) if c["start"]]
            hit = scan_hit(r)
            same += hit is not None
            other += len(calls) - (hit is not None)
            if hit and r["design"] == "microhomology":
                scan_mh += hit["microhomology"] == r["planted"]
        left = sum(1 for r in mh if int(r["planted_start"]) - int(r["start"]) == int(r["planted"]))
        mh_ok = sum(1 for r in mh if r["microhomology"] == r["planted"])
        false_rep = sum(1 for r in mh if r["flanking_repeat"] == "True")
        rep_ok = sum(1 for r in rp if r["flanking_repeat"] == "True")
        return (f"{label:<26} {same:>3}/{len(m3):<3}   {other:>3}     {scan_mh:>3}/{len(mh):<3}      {left:>3}/{len(mh):<3}"
                f"      {mh_ok:>3}/{len(mh):<3}       {false_rep:>3}/{len(mh):<3}       {rep_ok:>3}/{len(rp):<3}")

    print("Planted junctions, reads from a population at 30%")
    print(f"{'':<26} scan:     scan:   scan:        Module 3:    Module 3:     Module 3:     Module 3:")
    print(f"{'':<26} junction  other   microhom.    placed       microhom.     repeat where  planted repeat")
    print(f"{'':<26} found     calls   exact        leftmost     exact         none planted  called")
    for g in genomes:
        print(line(g, [r for r in rows if r["genome"] == g]))
    print(line("all", rows))
    m3 = [r for r in rows if r["source"] == "module3"]
    bad = [r for r in m3 if (
        (r["design"] == "microhomology" and (r["microhomology"] != r["planted"] or r["flanking_repeat"] == "True"))
        or (r["design"] == "repeat" and r["flanking_repeat"] != "True"))]
    if bad:
        print("\nModule 3 rows not as planted:")
        for r in bad:
            print(f"  {r['genome']} {r['id']}: planted {r['design']} {r['planted']} bp; reported repeat "
                  f"{r['flanking_repeat']} {r['repeat_len']} bp, microhomology {r['microhomology']}")
    rp = [r for r in m3 if r["design"] == "repeat" and r["flanking_repeat"] == "True"]
    if rp:
        print("\nReported length of the planted repeats (planted: reported):")
        for L in sorted({r["planted"] for r in rp}, key=int):
            print(f"  {L} bp: {sorted(int(r['repeat_len']) for r in rp if r['planted'] == L)}")

    # Position and microhomology at the planted repeats: the two copies of a repeat read
    # the same for some bases around the crossover, and the junction can lie anywhere
    # along them (expected_start and expected_microhomology, from the planted sequences).
    rp = [r for r in m3 if r["design"] == "repeat" and r.get("expected_start")]
    if rp:
        def as_expected(r, c):
            return c is not None and c["start"] == r["expected_start"] and \
                c["microhomology"] == r["expected_microhomology"]
        lengths = sorted({int(r["expected_microhomology"]) for r in rp})
        print(f"\nPlanted repeats: the junction can lie anywhere along the {'/'.join(map(str, lengths))} bp "
              "that the two copies share around the crossover.")
        print(f"  At the leftmost of these positions, with their number as microhomology: "
              f"Module 3 {sum(as_expected(r, r) for r in rp)}/{len(rp)}, "
              f"scan {sum(as_expected(r, scan_hit(r)) for r in rp)}/{len(rp)}")

    # Dosage: every junction was planted at the same share of cells.
    def dosage_line(label, rs):
        rs = [r for r in rs if r.get("dosage")]
        d = [float(r["dosage"]) for r in rs]
        if len(d) < 2:
            return None
        text = (f"  {label:<9} mean {statistics.mean(d):.3f}, sd {statistics.stdev(d):.3f}, "
                f"range {min(d):.3f}-{max(d):.3f}")
        ci = [r for r in rs if r.get("dosage_low") and r.get("dosage_high")]
        if ci:
            inside = sum(float(r["dosage_low"]) <= PLANTED_DOSAGE <= float(r["dosage_high"]) for r in ci)
            text += f"; 95% interval holds {PLANTED_DOSAGE:.2f} at {inside}/{len(ci)}"
        return text

    at_junction = [c for c in (scan_hit(r) for r in m3) if c is not None]
    lines = [dosage_line("Module 3:", m3), dosage_line("scan:", at_junction)]
    if any(lines):
        seeds = {r["read_seed"] for r in m3 if r.get("read_seed")}
        how = f"{len(m3)} read sets from {len(seeds)} seeds" if seeds else f"{len(m3)} read sets from one seed"
        print(f"\nDosage where {PLANTED_DOSAGE:.2f} was planted ({how})")
        for text in lines:
            if text:
                print(text)

    if rand:
        print("\nRandom junctions without repeats, described at exact coordinates (flanking repeat reported / junctions)")
        classes = []
        for r in rand:
            if r["size_class"] not in classes:
                classes.append(r["size_class"])
        print(f"{'':<26} " + "   ".join(f"{c:>11}" for c in classes))
        for g in sorted({r["genome"] for r in rand}) + ["all"]:
            cells = []
            for c in classes:
                rs = [r for r in rand if r["size_class"] == c and (g == "all" or r["genome"] == g)]
                cells.append(f"{sum(1 for r in rs if r['flanking_repeat'] == 'True'):>5}/{len(rs):<5}")
            print(f"{g:<26} " + "   ".join(cells))
        n = len(rand)
        mh = [int(r["microhomology"]) for r in rand]
        print(f"\nMicrohomology at the {n} random junctions: 0 bp {100 * mh.count(0) / n:.0f}%, "
              f"1 bp {100 * mh.count(1) / n:.0f}%, 2 bp {100 * mh.count(2) / n:.0f}%, "
              f"3 bp or more {100 * sum(1 for v in mh if v >= 3) / n:.0f}%, "
              f"5 bp or more {100 * sum(1 for v in mh if v >= 5) / n:.1f}%")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1])
