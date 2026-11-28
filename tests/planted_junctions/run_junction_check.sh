#!/bin/bash
# The junction that Modules 2 and 3 report, tested on planted junctions in each reference.
#
#   cd FOLDER_THAT_HOLDS genome/pop_ncbi_ref
#   THREADS=8 bash path/to/tests/planted_junctions/run_junction_check.sh
#
# For every <sp>/<sp>.fna under REF_DIR, and every file in EXTRA_REFS:
#   1. ten duplications with known junctions are planted (make_mutants.py);
#   2. for each, reads are simulated from a population carrying it at 30%, and the scan
#      (no coordinates) and Module 3 (coordinates given) are run on them (reads_planted.py);
#   3. the same description is tried at random junctions without repeats
#      (random_junctions.py, one background job per genome; RANDOM_N=0 skips it).
# Nothing in the reference folders is changed. Results go to OUT; OUT.tar.gz holds the
# tables and the summary, without genomes or reads.
#
# REF_DIR, EXTRA_REFS (space-separated FASTA files), SEED, RANDOM_N, OUT, THREADS and
# OVERWRITE can be set in the environment.

HERE=$(cd "$(dirname "$0")" && pwd)
THREADS=${THREADS:-8}
REF_DIR=${REF_DIR:-genome/pop_ncbi_ref}
EXTRA_REFS=${EXTRA_REFS:-}
SEED=${SEED:-1}
RANDOM_N=${RANDOM_N:-100}
VERSION=$(tandem --version | awk '{print $2}')
OUT=${OUT:-junction_check_v$VERSION}

if [ -e "$OUT" ] && [ "${OVERWRITE:-0}" != 1 ]; then
    echo "$OUT already exists; set OVERWRITE=1 to write into it" >&2
    exit 1
fi
mkdir -p "$OUT"

refs=()
for d in "$REF_DIR"/*/; do
    sp=$(basename "$d")
    [ -f "$d/$sp.fna" ] && refs+=("$d/$sp.fna")
done
for r in $EXTRA_REFS; do
    if [ -f "$r" ]; then refs+=("$r"); else echo "not found, skipped: $r"; fi
done
if [ ${#refs[@]} -eq 0 ]; then
    echo "no reference found under $REF_DIR or in EXTRA_REFS; run this from the folder that holds them" >&2
    exit 1
fi

echo "Tandem $VERSION, junction description on planted junctions: ${#refs[@]} genome(s) -> $OUT/"
if [ "$RANDOM_N" -gt 0 ]; then
    for ref in "${refs[@]}"; do
        sp=$(basename "$ref"); sp=${sp%.*}
        mkdir -p "$OUT/$sp"
        python "$HERE/random_junctions.py" "$ref" "$RANDOM_N" "$SEED" \
            > "$OUT/$sp/random_junctions.tsv" 2> "$OUT/$sp/random_junctions.log" &
    done
fi
for ref in "${refs[@]}"; do
    sp=$(basename "$ref"); sp=${sp%.*}
    if ! python "$HERE/make_mutants.py" "$ref" "$OUT/$sp/mutants" "$SEED"; then
        echo "FAILED: mutants for $sp" >&2
        continue
    fi
    python "$HERE/reads_planted.py" "$OUT/$sp" "$sp" "$THREADS" || echo "FAILED: reads test for $sp" >&2
done
if [ "$RANDOM_N" -gt 0 ]; then
    echo "waiting for the random-junction jobs"
    wait
fi

python "$HERE/summarise_reads.py" "$OUT" | tee "$OUT/summary.txt"
tar czf "$OUT.tar.gz" -C "$(dirname "$OUT")" --exclude='*.fna' --exclude='*.fq' --exclude='work' \
    "$(basename "$OUT")"
echo "Results: $OUT.tar.gz"
