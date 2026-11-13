"""Module 2 scan checks on simulated reads with known answers.

    python tests/check_scan.py WORKDIR

Covers the cases real data throws at the scanner: reads running across the
origin of a circular chromosome, a duplication on a second contig, a small
indel that should stay below the size floor, a repeat right beside the
junction, and background marked from a control sample. Prints PASS/FAIL per
case and exits non-zero if anything failed.
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

    failed = 0
    for name, ok, detail in results:
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<42} {detail}")
        failed += not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "/tmp/tandem_scan_checks"))
