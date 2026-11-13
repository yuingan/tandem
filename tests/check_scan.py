"""Module 2 scan checks on simulated reads with known answers.

    python tests/check_scan.py WORKDIR

Covers the cases real data throws at the scanner: reads running across the
origin of a circular chromosome, a duplication on a second contig, a small
indel that should stay below the size floor, a repeat right beside the
junction, background marked from a control sample, boundaries given only
roughly, and a junction no read can cross. Prints PASS/FAIL per case and
exits non-zero if anything failed.
"""

import csv
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SIMULATE = os.path.join(HERE, "simulate_reads.py")


def write_fasta(path, records):
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                f.write(seq[i:i + 70] + "\n")


def simulate(ref, prefix, start, end, freq, seed, depth=100):
    subprocess.run(
        [sys.executable, SIMULATE, ref, prefix, "--start", str(start),
         "--end", str(end), "--frequency", str(freq), "--depth", str(depth),
         "--circular", "--seed", str(seed)],
        check=True, capture_output=True)
    return prefix + "_R1.fq", prefix + "_R2.fq"


def concat(paths, out):
    with open(out, "w") as o:
        for p in paths:
            with open(p) as f:
                o.write(f.read())
    return out


def run_scan(ref, r1, r2, outdir, extra=()):
    subprocess.run(["tandem", "-r", ref, "-i", r1, "-I", r2, "-iso",
                    "-o", outdir, *extra], check=True, capture_output=True)
    with open(os.path.join(outdir, "junctions.tsv")) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def run_population(ref, r1, r2, metadata, outdir):
    subprocess.run(["tandem", "-r", ref, "-i", r1, "-I", r2, "-pop", "-m", metadata,
                    "-o", outdir], check=True, capture_output=True)
    with open(os.path.join(outdir, "junction_quantification.tsv")) as f:
        lines = [line for line in f if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def table_text(outdir):
    """The summary table of the HTML report, tags stripped."""
    import re
    with open(os.path.join(outdir, "tandem_report.html")) as f:
        html = f.read()
    body = html.split("<tbody>", 1)[1].split("</tbody>", 1)[0] if "<tbody>" in html else ""
    return " ".join(re.sub(r"<[^>]+>", " ", body).split())


def near(row, seq_id, start, end, slop=5):
    return (row["seq_id"] == seq_id
            and abs(int(row["start"]) - start) <= slop
            and abs(int(row["end"]) - end) <= slop)


def main(work):
    os.makedirs(work, exist_ok=True)
    random.seed(5)
    chrom = "".join(random.choice("ACGT") for _ in range(200_000))
    plasmid = "".join(random.choice("ACGT") for _ in range(80_000))

    chrom_fa = os.path.join(work, "chrom.fna")
    plasmid_fa = os.path.join(work, "plasmid.fna")
    both_fa = os.path.join(work, "chrom_plasmid.fna")
    write_fasta(chrom_fa, [("chrom", chrom)])
    write_fasta(plasmid_fa, [("plasmid", plasmid)])
    write_fasta(both_fa, [("chrom", chrom), ("plasmid", plasmid)])

    results = []

    # 1. Circular chromosome: the 10 kb duplication and nothing else.
    r1, r2 = simulate(chrom_fa, os.path.join(work, "circ"), 50000, 60000, 0.4, 11)
    rows = run_scan(chrom_fa, r1, r2, os.path.join(work, "out_circ"))
    ok = len(rows) == 1 and near(rows[0], "chrom", 50000, 60000)
    dosage = rows[0]["dosage"] if rows else "-"
    results.append(("circular chromosome, 10 kb at f=0.4", ok,
                    f"{len(rows)} junction(s), dosage {dosage}"))
    circ = rows

    # The report shows what junctions.tsv holds, not the old module 2 fields.
    text = table_text(os.path.join(work, "out_circ"))
    ok = "50,000" in text and "60,000" in text and (not circ or circ[0]["junction_reads"] in text.split())
    results.append(("isolate report lists the junction", ok, text[:60]))

    # 2. Two contigs, duplication on the second one.
    c1, c2 = simulate(chrom_fa, os.path.join(work, "wt"), 1, 2, 0.0, 12)
    p1, p2 = simulate(plasmid_fa, os.path.join(work, "pdup"), 30000, 40000, 0.5, 13)
    m1 = concat([c1, p1], os.path.join(work, "mix_R1.fq"))
    m2 = concat([c2, p2], os.path.join(work, "mix_R2.fq"))
    rows = run_scan(both_fa, m1, m2, os.path.join(work, "out_two"))
    ok = len(rows) == 1 and near(rows[0], "plasmid", 30000, 40000)
    dosage = rows[0]["dosage"] if rows else "-"
    results.append(("second contig, 10 kb at f=0.5", ok,
                    f"{len(rows)} junction(s), dosage {dosage}"))

    # 3. A 3 bp duplication is an indel, not something to report.
    s1, s2 = simulate(chrom_fa, os.path.join(work, "small"), 80000, 80002, 1.0, 14)
    rows = run_scan(chrom_fa, s1, s2, os.path.join(work, "out_small"))
    results.append(("3 bp duplication stays below --min-size", len(rows) == 0,
                    f"{len(rows)} junction(s)"))

    # 4. The first 80 bp of the unit also occur elsewhere, so most reads across
    #    the junction can't be placed. The raw ratio drops well below the
    #    truth; the corrected one should not.
    random.seed(12)
    base = list("".join(random.choice("ACGT") for _ in range(300_000)))
    base[200000:200080] = base[100000:100080]
    rep_fa = os.path.join(work, "repeat_at_start.fna")
    write_fasta(rep_fa, [("chrom", "".join(base))])
    # deep enough that the reads left countable (about 40%) pin the estimate
    r1, r2 = simulate(rep_fa, os.path.join(work, "rep"), 100001, 110000, 0.5, 15,
                      depth=600)
    rows = [r for r in run_scan(rep_fa, r1, r2, os.path.join(work, "out_rep"))
            if near(r, "chrom", 100001, 110000)]
    ok = (len(rows) == 1 and float(rows[0]["dosage_raw"]) < 0.3
          and 0.3 <= float(rows[0]["dosage"]) <= 0.7)
    detail = (f"raw {rows[0]['dosage_raw']}, corrected {rows[0]['dosage']} "
              f"(truth 0.5)") if rows else "not found"
    results.append(("repeat beside the junction, f=0.5", ok, detail))

    # 5. A junction that is also in the control is marked with its support there.
    c1, c2 = simulate(chrom_fa, os.path.join(work, "ctrl"), 50000, 60000, 0.2, 16)
    ctrl_rows = run_scan(chrom_fa, c1, c2, os.path.join(work, "out_ctrl"),
                         extra=["--min-support", "1"])
    rows = run_scan(chrom_fa, os.path.join(work, "circ_R1.fq"),
                    os.path.join(work, "circ_R2.fq"), os.path.join(work, "out_vs_ctrl"),
                    extra=["--control-junctions",
                           os.path.join(work, "out_ctrl", "junctions.tsv")])
    ok = (len(ctrl_rows) >= 1 and len(rows) == 1
          and int(rows[0]["control_reads"]) > 0)
    detail = f"control_reads {rows[0]['control_reads']}" if rows else "not found"
    results.append(("junction also in the control is marked", ok, detail))

    # 6. Population mode at known coordinates, on the reads from case 1: the
    #    junction as given, the same 6 bp off, and a stretch with no
    #    duplication. The first two are one junction and must agree.
    metadata = os.path.join(work, "known.tsv")
    with open(metadata, "w") as f:
        f.write("given\t50000\t60000\nshifted\t50006\t60006\nabsent\t120000\t125000\n")
    rows = {r["name"]: r for r in run_population(
        chrom_fa, os.path.join(work, "circ_R1.fq"), os.path.join(work, "circ_R2.fq"),
        metadata, os.path.join(work, "out_pop"))}
    given, shifted, absent = rows.get("given"), rows.get("shifted"), rows.get("absent")
    ok = (given is not None and shifted is not None and absent is not None
          and int(given["junction_reads"]) > 0
          and 0.2 <= float(given["dosage"]) <= 0.6
          and (shifted["dup_start"], shifted["dosage"]) == (given["dup_start"], given["dosage"])
          and int(absent["junction_reads"]) == 0 and float(absent["dosage"]) == 0)
    detail = (f"dosage {given['dosage']} (truth 0.4), off-by-6 {shifted['dosage']}, "
              f"absent {absent['dosage']}") if ok or (given and shifted and absent) else "missing rows"
    results.append(("population mode at known junctions", ok, detail))

    # 7. Isolate mode with boundaries read off a coverage plot, 150 bp out: every
    #    start paired with every end after it, each reported once.
    out = os.path.join(work, "out_iso_given")
    subprocess.run(["tandem", "-r", chrom_fa, "-i", os.path.join(work, "circ_R1.fq"),
                    "-I", os.path.join(work, "circ_R2.fq"), "-iso",
                    "-s", "49850", "120000", "-e", "60150", "125000", "-o", out],
                   check=True, capture_output=True)
    with open(os.path.join(out, "junctions.tsv")) as f:
        rows = {(int(r["provided_start"]), int(r["provided_end"])): r
                for r in csv.DictReader(f, delimiter="\t")}
    hit = rows.get((49850, 60150))
    others = [r for key, r in rows.items() if key != (49850, 60150)]
    ok = (hit is not None and near(hit, "chrom", 50000, 60000)
          and circ and hit["junction_reads"] == circ[0]["junction_reads"]
          and len(others) == 2 and all(int(r["junction_reads"]) == 0 for r in others))
    detail = (f"{hit['start']}-{hit['end']}, {hit['junction_reads']} reads; "
              f"{len(others)} other pair(s) without reads") if hit else f"{len(rows)} row(s)"
    results.append(("isolate mode with -s/-e, 150 bp out", ok, detail))

    # 8. The first 400 bp of the unit also occur elsewhere, longer than a read:
    #    no read can be placed across the junction, so the dosage comes from
    #    depth.
    random.seed(13)
    base = list("".join(random.choice("ACGT") for _ in range(300_000)))
    base[200000:200400] = base[100000:100400]
    hidden_fa = os.path.join(work, "repeat_400_at_start.fna")
    write_fasta(hidden_fa, [("chrom", "".join(base))])
    h1, h2 = simulate(hidden_fa, os.path.join(work, "hidden"), 100001, 110000, 0.5, 17)
    metadata = os.path.join(work, "hidden.tsv")
    with open(metadata, "w") as f:
        f.write("hidden\t100001\t110000\n")
    rows = run_population(hidden_fa, h1, h2, metadata, os.path.join(work, "out_hidden"))
    r = rows[0] if rows else None
    ok = (r is not None and int(r["junction_reads"]) == 0 and float(r["junction_visible"]) == 0
          and r["dosage_from"] == "depth" and 0.35 <= float(r["dosage"]) <= 0.65)
    detail = (f"junction visible {r['junction_visible']}, dosage {r['dosage']} from "
              f"{r['dosage_from']} (truth 0.5)") if r else "no row"
    results.append(("junction inside a repeat: depth instead", ok, detail))
    hidden = r

    # 9. Every estimate comes with an interval around it: from the read counts
    #    for a junction (from 0 where no read crosses it), from the noise of
    #    depth where the dosage comes from depth.
    def spans(row):
        return (row is not None and row["dosage_low"] != "" and row["dosage_high"] != ""
                and float(row["dosage_low"]) <= float(row["dosage"]) <= float(row["dosage_high"]))
    ok = (spans(given) and spans(absent) and float(absent["dosage_low"]) == 0
          and float(absent["dosage_high"]) > 0 and spans(hidden)
          and float(hidden["depth_noise"]) > 0)
    detail = (f"{given['dosage_low']}-{given['dosage_high']} around {given['dosage']}; "
              f"depth noise {hidden['depth_noise']}") if given and hidden else "missing rows"
    results.append(("intervals around the estimates", ok, detail))

    failed = 0
    for name, ok, detail in results:
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<42} {detail}")
        failed += not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "/tmp/tandem_scan_checks"))
