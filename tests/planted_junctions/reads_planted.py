#!/usr/bin/env python3
"""The junction description from reads, tested on planted junctions.

For each duplication planted by make_mutants.py, the reference state (one copy of the
unit) is cut 120 kb either side of it, reads are simulated from a population in which
30% of cells carry the duplication (150 bp pairs, 100x, 0.1% errors), and Tandem is run
twice on them: the scan, with no coordinates, and Module 3, with the coordinates given.
The junction each reports and its flanking-repeat and microhomology columns are written
next to what was planted.

Every junction gets its own read seed, made from the genome's label and the junction
(read_seed in the table). The simulator's draws do not depend on the sequence, so with
one seed for all, two excerpts of the same length receive their fragments at the same
positions, and units of similar size have their junction in the same stretch of the same
layout: the read counts of different genomes would then repeat each other.

Usage: reads_planted.py GENOME_DIR LABEL [THREADS]
  GENOME_DIR holds mutants/truth.tsv and mutants/mNN.fna; writes GENOME_DIR/reads_planted.tsv.
Needs tandem on the PATH and the suite's read simulator, tests/simulate_reads.py (one folder up).
"""
import csv
import os
import shutil
import subprocess
import sys
import zlib

FLANK = 120_000
HERE = os.path.dirname(os.path.abspath(__file__))
SIMULATE = next((p for p in (os.path.join(HERE, "simulate_reads.py"),
                             os.path.join(os.path.dirname(HERE), "simulate_reads.py"))
                 if os.path.exists(p)), None)
COLUMNS = ["genome", "id", "design", "planted", "unit", "planted_start", "planted_end", "source",
           "start", "end", "junction_reads", "reference_reads", "dosage", "flanking_repeat",
           "repeat_len", "repeat_identity", "microhomology", "dosage_low", "dosage_high",
           "expected_start", "expected_microhomology", "read_seed"]


def read_record(path, name):
    recs, cur, seq = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    recs[cur] = "".join(seq).upper()
                cur, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if cur is not None:
        recs[cur] = "".join(seq).upper()
    return recs[name] if name in recs else max(recs.values(), key=len)


def run(cmd):
    done = subprocess.run(cmd, capture_output=True, text=True)
    if done.returncode != 0:
        sys.exit(f"FAILED: {' '.join(cmd)}\n{done.stderr[-2000:]}")


def main(base, label, threads="4"):
    work = os.path.join(base, "work")
    os.makedirs(work, exist_ok=True)
    out = os.path.join(base, "reads_planted.tsv")
    with open(out, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")

        def write(row):
            fh.write("\t".join(str(row.get(c, "")) for c in COLUMNS) + "\n")

        for t in csv.DictReader(open(os.path.join(base, "mutants", "truth.tsv")), delimiter="\t"):
            M = read_record(os.path.join(base, "mutants", f"{t['id']}.fna"), t.get("seq_id"))
            c1s, c1e, c2e = int(t["copy1_start"]), int(t["copy1_end"]), int(t["copy2_end"])
            ref = M[:c1e] + M[c2e:]                  # reference state: the second copy removed
            s0, e0 = c1s - 1, c1e
            if t["design"] == "repeat":
                n = int(t["repeat_len"])
                x = n // 2
                s0, e0 = s0 + x, e0 + x              # the unit at the planted crossover
                design, planted = "repeat", t["repeat_len"]
                # The two copies of the repeat read the same for some bases either side of
                # the crossover, and the junction can lie anywhere along them: it is
                # expected at the leftmost of these positions, with their number as its
                # microhomology.
                R1, R2 = M[c1s - 1:c1s - 1 + n], M[c2e:c2e + n]
                before = next((i for i in range(x) if R1[x - 1 - i] != R2[x - 1 - i]), x)
                after = next((i for i in range(n - x) if R1[x + i] != R2[x + i]), n - x)
                shift, expected_mh = before, before + after
            else:
                design, planted = "microhomology", t["microhomology"]
                shift = expected_mh = int(planted)
            lo, hi = max(0, s0 - FLANK), min(len(ref), e0 + FLANK)
            ex = ref[lo:hi]
            a, b = s0 - lo + 1, e0 - lo              # the unit in the excerpt, 1-based inclusive
            fa = os.path.join(work, "ref.fna")
            with open(fa, "w") as o:
                o.write(">excerpt\n" + "\n".join(ex[i:i + 70] for i in range(0, len(ex), 70)) + "\n")
            seed = zlib.crc32(f"{label} {t['id']} {c1s} {c1e}".encode())
            prefix = os.path.join(work, "reads")
            run([sys.executable, SIMULATE, fa, prefix, "--start", str(a),
                 "--end", str(b), "--frequency", "0.3", "--depth", "100", "--seed", str(seed)])
            r1, r2 = prefix + "_R1.fq", prefix + "_R2.fq"
            meta = os.path.join(work, "known.tsv")
            with open(meta, "w") as o:
                o.write(f"{t['id']}\t{a}\t{b}\n")
            head = {"genome": label, "id": t["id"], "design": design, "planted": planted,
                    "unit": b - a + 1, "planted_start": a, "planted_end": b,
                    "expected_start": a - shift, "expected_microhomology": expected_mh,
                    "read_seed": seed}

            od = os.path.join(work, "pop")           # Module 3, coordinates given
            shutil.rmtree(od, ignore_errors=True)
            run(["tandem", "-r", fa, "-i", r1, "-I", r2, "-pop", "-m", meta, "-o", od, "-t", threads,
                 "--no-circular"])
            rows = [l for l in open(os.path.join(od, "junction_quantification.tsv")) if not l.startswith("#")]
            for r in csv.DictReader(rows, delimiter="\t"):
                write(dict(head, source="module3", start=r["dup_start"], end=r["dup_end"],
                           junction_reads=r["junction_reads"], reference_reads=r["reference_reads"],
                           dosage=r["dosage"], dosage_low=r["dosage_low"], dosage_high=r["dosage_high"],
                           flanking_repeat=r["flanking_repeat"], repeat_len=r["repeat_len"],
                           repeat_identity=r["repeat_identity"], microhomology=r["microhomology_bp"]))

            od = os.path.join(work, "scan")          # the scan, no coordinates
            shutil.rmtree(od, ignore_errors=True)
            run(["tandem", "-r", fa, "-i", r1, "-I", r2, "-iso", "-o", od, "-t", threads, "--no-circular"])
            found = list(csv.DictReader(open(os.path.join(od, "junctions.tsv")), delimiter="\t"))
            if not found:
                write(dict(head, source="scan", junction_reads=0))
            for r in found:
                write(dict(head, source="scan", start=r["start"], end=r["end"],
                           junction_reads=r["junction_reads"], reference_reads=r["reference_reads"],
                           dosage=r["dosage"], dosage_low=r["dosage_low"], dosage_high=r["dosage_high"],
                           microhomology=r["microhomology"]))
            fh.flush()
    shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    if SIMULATE is None:
        sys.exit("simulate_reads.py not found beside this script or one folder up")
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "4")
