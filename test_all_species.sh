#!/bin/bash
# Test Tandem module 1 on all reference genomes found in genome/pop_ncbi_ref/
# Run from the working directory containing genome/pop_ncbi_ref/

set -euo pipefail

REF_DIR="genome/pop_ncbi_ref"
OUT_BASE="test_module1_all_species"
THREADS=8

mkdir -p "$OUT_BASE"

# Optional human-readable names. Any species not listed here falls back to
# its directory name. Add the 7th species' pretty label below if you want one.
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

# Helper: return the display name (pretty label, or directory name as fallback)
disp() { echo "${SPECIES[$1]:-$1}"; }

echo "============================================"
echo "Tandem Module 1: Testing ${#SPECIES_LIST[@]} species"
echo "============================================"
echo ""

# Run each species
for sp in "${SPECIES_LIST[@]}"; do
    REF="$REF_DIR/$sp/${sp}.fna"
    OUTDIR="$OUT_BASE/$sp"

    echo "--------------------------------------------"
    echo "Processing: $(disp "$sp") ($sp)"
    echo "--------------------------------------------"

    tandem -r "$REF" -o "$OUTDIR" -t "$THREADS" -v

    echo ""
done

# Summary across all species
echo ""
echo "============================================"
echo "Summary across all species"
echo "============================================"
echo ""
printf "%-30s %8s %8s %8s %8s %8s\n" "Species" "Tandems" "HR" "non-HR" "Flagged" "Circ"
printf "%-30s %8s %8s %8s %8s %8s\n" "-------" "-------" "----" "------" "-------" "----"

for sp in "${SPECIES_LIST[@]}"; do
    TSV="$OUT_BASE/$sp/tandem_duplications.tsv"
    if [ -f "$TSV" ]; then
        total=$(tail -n +2 "$TSV" | wc -l)
        hr=$(tail -n +2 "$TSV" | awk -F'\t' '$14=="True"' | wc -l)
        nonhr=$((total - hr))
        flagged=$(tail -n +2 "$TSV" | awk -F'\t' '$22!=""' | wc -l)
        circular=$(tail -n +2 "$TSV" | awk -F'\t' '$8=="True"' | wc -l)
        printf "%-30s %8d %8d %8d %8d %8d\n" "$(disp "$sp")" "$total" "$hr" "$nonhr" "$flagged" "$circular"
    else
        printf "%-30s %8s %8s %8s %8s %8s\n" "$(disp "$sp")" "N/A" "N/A" "N/A" "N/A" "N/A"
    fi
done

echo ""
echo "Microhomology distribution (non-HR, all species combined):"
echo ""
printf "%6s %8s\n" "MH(bp)" "Count"
printf "%6s %8s\n" "------" "-----"

# Combine all TSVs and count microhomology for non-HR
for sp in "${SPECIES_LIST[@]}"; do
    TSV="$OUT_BASE/$sp/tandem_duplications.tsv"
    if [ -f "$TSV" ]; then
        tail -n +2 "$TSV" | awk -F'\t' '$14=="False"'
    fi
done | awk -F'\t' '{print $18}' | sort -n | uniq -c | sort -rn | head -20 | while read count mh; do
    printf "%6s %8s\n" "$mh" "$count"
done

echo ""
echo "Results in: $OUT_BASE/"
echo "HTML reports: $OUT_BASE/*/tandem_report.html"
