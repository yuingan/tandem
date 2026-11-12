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
# Ask the installed command, not `python -c "import tandem"`: run from a
# folder holding a tandem/ source tree, the import picks up that copy instead.
VERSION=$(tandem --version | awk '{print $2}')
OUT_BASE=${OUT_BASE:-test_module1_v$VERSION}
THREADS=${THREADS:-8}

# Never write over an earlier run unless asked to.
if [ -n "$(ls -A "$OUT_BASE" 2>/dev/null)" ] && [ "${OVERWRITE:-0}" != 1 ]; then
    echo "ERROR: $OUT_BASE already holds results. Set OUT_BASE to a new folder," >&2
    echo "       or OVERWRITE=1 to replace them." >&2
    exit 1
fi

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

# One genome failing shouldn't cost the others; failures are listed at the end.
FAILED=()
for sp in "${SPECIES_LIST[@]}"; do
    REF="$REF_DIR/$sp/${sp}.fna"
    OUTDIR="$OUT_BASE/$sp"

    echo "--------------------------------------------"
    echo "Processing: $(disp "$sp") ($sp)"
    echo "--------------------------------------------"

    if ! tandem -r "$REF" -o "$OUTDIR" -t "$THREADS" -v; then
        echo "FAILED: $sp" >&2
        FAILED+=("$sp")
    fi

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

if [ "${#FAILED[@]}" -gt 0 ]; then
    echo ""
    echo "FAILED: ${FAILED[*]}" >&2
    exit 1
fi
