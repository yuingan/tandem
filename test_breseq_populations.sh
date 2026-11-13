#!/bin/bash
# breseq on the artificial populations, scored the same way as Module 3.
#
#   cd /mnt/d/tandem_software
#   TANDEM=module3_v0.4.0.dev13_benchmark JOBS=5 THREADS=8 \
#       bash ~/tandem_git/test_breseq_populations.sh
#
# MODE=given (the default) passes breseq every duplication's junction with
# --user-evidence-gd, as Module 3 measured it in TANDEM, so breseq reports
# support for each whatever its frequency. MODE=denovo runs breseq -p alone,
# to see which it finds unprompted. SPECIES and DAYS narrow the run, e.g.
# SPECIES=ecoli_K12_MG1655 DAYS=40 for a first try.
#
# breseq's junction frequency is J / (J + S): for a duplication carried by a
# share p of cells that is p/(1+p), so it is converted back as F/(1-F).
# JOBS runs of THREADS each go at once. TIDY=0 keeps breseq's intermediate
# folders, which otherwise go once a run has been converted.

set -euo pipefail
shopt -s nullglob

HERE=$(cd "$(dirname "$0")" && pwd)
ART_DIR=${ART_DIR:-genome/artificial_pop_seqs}
REF_DIR=${REF_DIR:-genome/pop_ncbi_ref}
MODE=${MODE:-given}
JOBS=${JOBS:-4}
THREADS=${THREADS:-8}
TIDY=${TIDY:-1}
TANDEM=${TANDEM:?set TANDEM to a Module 3 benchmark folder, e.g. module3_v0.4.0.dev13_benchmark}
VERSION=$(breseq --version 2>&1 | awk 'NR==1 {print $2}')
OUT=${OUT:-breseq_v${VERSION}_${MODE}_benchmark}

case $MODE in given|denovo) ;; *) echo "MODE must be given or denovo" >&2; exit 1 ;; esac
mkdir -p "$OUT/logs" "$OUT/evidence" "$OUT/runs"
echo "breseq $VERSION, $MODE, on $ART_DIR -> $OUT/"

tasks=()
for gt in "$ART_DIR"/*/ground_truth.json; do
    species=$(basename "$(dirname "$gt")")
    [ -n "${SPECIES:-}" ] && [[ " $SPECIES " != *" $species "* ]] && continue
    python "$HERE/tests/breseq_junctions.py" evidence --tandem-dir "$TANDEM" \
        --species "$species" --gd "$OUT/evidence/$species.gd" \
        --targets "$OUT/evidence/${species}_targets.tsv"
    for day_dir in "$ART_DIR/$species"/day_*; do
        day=${day_dir##*_}
        [ -n "${DAYS:-}" ] && [[ " $DAYS " != *" $day "* ]] && continue
        tasks+=("$species $day")
    done
done
if [ ${#tasks[@]} -eq 0 ]; then
    echo "Nothing to run: check SPECIES, DAYS and ART_DIR" >&2
    exit 1
fi

run_one() {
    local species=$1 day=$2
    local refs=("$REF_DIR/$species"*/*.fna)
    local reads=("$ART_DIR/$species/day_$day"/*_R[12].fastq*)
    local name=${species}_day${day}
    local run=$OUT/runs/$name
    local extra=()
    [ "$MODE" = given ] && extra=(--user-evidence-gd "$OUT/evidence/$species.gd")
    SECONDS=0
    if breseq -p -j "$THREADS" -r "${refs[0]}" -o "$run" "${extra[@]}" "${reads[@]}" \
            > "$OUT/logs/$name.log" 2>&1 \
       && python "$HERE/tests/breseq_junctions.py" convert --breseq-dir "$run" \
            --targets "$OUT/evidence/${species}_targets.tsv" --mode "$MODE" \
            --out "$OUT/$name/junction_quantification.tsv" >> "$OUT/logs/$name.log" 2>&1; then
        echo "  $name: done in ${SECONDS}s"
        [ "$TIDY" = 1 ] && rm -rf "$run"/0[0-9]_*
    else
        echo "  $name: FAILED, see $OUT/logs/$name.log" >&2
        touch "$OUT/logs/$name.failed"
    fi
    return 0
}
export -f run_one
export ART_DIR REF_DIR MODE THREADS TIDY OUT HERE

rm -f "$OUT"/logs/*.failed
printf '%s\n' "${tasks[@]}" | xargs -P "$JOBS" -L 1 bash -c 'run_one $0 $1'

echo
python "$HERE/tests/score_artificial_populations.py" --gt-dir "$ART_DIR" \
    --results-dir "$OUT" | tee "$OUT/summary.txt"

failed=("$OUT"/logs/*.failed)
if [ ${#failed[@]} -gt 0 ]; then
    echo "FAILED: ${failed[*]}" >&2
    exit 1
fi
