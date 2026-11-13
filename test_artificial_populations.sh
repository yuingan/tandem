#!/bin/bash
# Module 3 on the artificial populations behind Figure 3, then the scores.
#
#   cd /mnt/d/tandem_software
#   THREADS=44 bash ~/tandem_git/test_artificial_populations.sh
#
# Every species folder under ART_DIR with a ground_truth.json is run, one
# tandem -pop per day_<d> folder, which holds the reads and metadata.tsv.
# The reference is the .fna under REF_DIR/<species>*/. Results go to
# module3_v<version>_benchmark/, named after the installed tandem so an old
# run is never overwritten.

set -euo pipefail
shopt -s nullglob

HERE=$(cd "$(dirname "$0")" && pwd)
THREADS=${THREADS:-8}
ART_DIR=${ART_DIR:-genome/artificial_pop_seqs}
REF_DIR=${REF_DIR:-genome/pop_ncbi_ref}
VERSION=$(tandem --version | awk '{print $2}')
OUT=${OUT:-module3_v${VERSION}_benchmark}

if [ -e "$OUT" ] && [ "${OVERWRITE:-0}" != 1 ]; then
    echo "$OUT already exists; set OVERWRITE=1 to write into it" >&2
    exit 1
fi
mkdir -p "$OUT/logs"

echo "Tandem $VERSION, Module 3 on $ART_DIR -> $OUT/"
FAILED=()
for gt in "$ART_DIR"/*/ground_truth.json; do
    species=$(basename "$(dirname "$gt")")
    refs=("$REF_DIR/$species"*/*.fna)
    if [ ${#refs[@]} -eq 0 ]; then
        echo "SKIP $species: no .fna under $REF_DIR/$species*/"
        continue
    fi
    ref=${refs[0]}
    echo "$species  (reference $ref)"

    for day_dir in "$ART_DIR/$species"/day_*; do
        day=${day_dir##*_}
        r1=("$day_dir"/*_R1.fastq*)
        r2=("$day_dir"/*_R2.fastq*)
        meta="$day_dir/metadata.tsv"
        if [ ${#r1[@]} -eq 0 ] || [ ! -f "$meta" ]; then
            echo "  day $day: SKIP, reads or metadata.tsv missing"
            continue
        fi
        pair=()
        [ ${#r2[@]} -gt 0 ] && pair=(-I "${r2[0]}")
        log="$OUT/logs/${species}_day${day}.log"
        SECONDS=0
        if tandem -r "$ref" -i "${r1[0]}" "${pair[@]}" -pop -m "$meta" -t "$THREADS" \
                -o "$OUT/${species}_day${day}" > "$log" 2>&1; then
            echo "  day $day: done in ${SECONDS}s"
        else
            echo "  day $day: FAILED, see $log" >&2
            FAILED+=("${species}_day${day}")
        fi
    done
done

echo
python "$HERE/tests/score_artificial_populations.py" --gt-dir "$ART_DIR" \
    --results-dir "$OUT" | tee "$OUT/summary.txt"

if [ ${#FAILED[@]} -gt 0 ]; then
    echo "FAILED: ${FAILED[*]}" >&2
    exit 1
fi
