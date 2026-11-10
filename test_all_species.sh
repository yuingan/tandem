#!/bin/bash
# Run Tandem module 1 on every reference genome under genome/pop_ncbi_ref/
# and summarise. Run from the folder that contains genome/pop_ncbi_ref/.
#
#   bash path/to/test_all_species.sh
#
# REF_DIR, OUT_BASE and THREADS can be set in the environment. OUT_BASE
# defaults to a folder named after the installed Tandem version, so running a
# new version never overwrites results from an older one.

set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
REF_DIR=${REF_DIR:-genome/pop_ncbi_ref}
VERSION=$(python -c "import tandem; print(tandem.__version__)")
OUT_BASE=${OUT_BASE:-test_module1_v$VERSION}
THREADS=${THREADS:-8}

mkdir -p "$OUT_BASE"

# Optional human-readable names. Any species not listed here falls back to
# its directory name.
declare -A SPECIES=(
    ["ecoli_K12_MG1655"]="E. coli K-12 MG1655"
    ["abaumannii_ATCC17978"]="A. baumannii ATCC 17978"
    ["scoelicolor_A3_2"]="S. coelicolor A3(2)"
    ["mtuberculosis_H37Rv"]="M. tuberculosis H37Rv"
    ["mpneumoniae_M129"]="M. pneumoniae M129"
    ["paeruginosa_PAO1"]="P. aeruginosa PAO1"
    ["mgenitalium_G37"]="M. genitalium G37"
)

# Single source of truth: auto-discover every <sp>/<sp>.fna under REF_DIR.
# Sorted for stable, reproducible ordering.
SPECIES_LIST=()
for d in "$REF_DIR"/*/; do
    sp=$(basename "$d")
    if [ -f "$d/${sp}.fna" ]; then
        SPECIES_LIST+=("$sp")
    fi
done
IFS=$'\n' SPECIES_LIST=($(sort <<<"${SPECIES_LIST[*]}")); unset IFS

if [ "${#SPECIES_LIST[@]}" -eq 0 ]; then
    echo "ERROR: no <sp>/<sp>.fna references found under $REF_DIR" >&2
    exit 1
fi

disp() { echo "${SPECIES[$1]:-$1}"; }

echo "============================================"
echo "Tandem $VERSION module 1: ${#SPECIES_LIST[@]} species -> $OUT_BASE/"
echo "============================================"
echo ""

for sp in "${SPECIES_LIST[@]}"; do
    REF="$REF_DIR/$sp/${sp}.fna"
    OUTDIR="$OUT_BASE/$sp"

    echo "--------------------------------------------"
    echo "Processing: $(disp "$sp") ($sp)"
    echo "--------------------------------------------"

    tandem -r "$REF" -o "$OUTDIR" -t "$THREADS" -v

    echo ""
done

echo "============================================"
echo "Summary"
echo "============================================"
python "$HERE/tests/compare_module1.py" "$OUT_BASE"

echo ""
echo "Results in: $OUT_BASE/"
echo "HTML reports: $OUT_BASE/*/tandem_report.html"
echo "Compare with an earlier run:"
echo "  python $HERE/tests/compare_module1.py OLD_RUN $OUT_BASE"
