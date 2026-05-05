#!/bin/bash
# Test Tandem module 1 on all 6 reference genomes
# Run from the working directory containing genome/pop_ncbi_ref/

set -euo pipefail

REF_DIR="genome/pop_ncbi_ref"
OUT_BASE="test_module1_all_species"
THREADS=8

mkdir -p "$OUT_BASE"

# Define species
declare -A SPECIES=(
    ["ecoli_K12_MG1655"]="E. coli K-12 MG1655"
    ["abaumannii_ATCC17978"]="A. baumannii ATCC 17978"
    ["scoelicolor_A3_2"]="S. coelicolor A3(2)"
    ["mtuberculosis_H37Rv"]="M. tuberculosis H37Rv"
    ["mpneumoniae_M129"]="M. pneumoniae M129"
    ["paeruginosa_PAO1"]="P. aeruginosa PAO1"
)

echo "============================================"
echo "Tandem Module 1: Testing all 6 species"
echo "============================================"
echo ""

# Run each species
for sp in ecoli_K12_MG1655 abaumannii_ATCC17978 scoelicolor_A3_2 mtuberculosis_H37Rv mpneumoniae_M129 paeruginosa_PAO1; do
    REF="$REF_DIR/$sp/${sp}.fna"
    OUTDIR="$OUT_BASE/$sp"

    if [ ! -f "$REF" ]; then
        echo "[SKIP] $sp: reference not found at $REF"
        continue
    fi

    echo "--------------------------------------------"
    echo "Processing: ${SPECIES[$sp]} ($sp)"
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

for sp in ecoli_K12_MG1655 abaumannii_ATCC17978 scoelicolor_A3_2 mtuberculosis_H37Rv mpneumoniae_M129 paeruginosa_PAO1; do
    TSV="$OUT_BASE/$sp/tandem_duplications.tsv"
    if [ -f "$TSV" ]; then
        total=$(tail -n +2 "$TSV" | wc -l)
        hr=$(tail -n +2 "$TSV" | awk -F'\t' '$14=="True"' | wc -l)
        nonhr=$((total - hr))
        flagged=$(tail -n +2 "$TSV" | awk -F'\t' '$22!=""' | wc -l)
        circular=$(tail -n +2 "$TSV" | awk -F'\t' '$8=="True"' | wc -l)
        printf "%-30s %8d %8d %8d %8d %8d\n" "${SPECIES[$sp]}" "$total" "$hr" "$nonhr" "$flagged" "$circular"
    else
        printf "%-30s %8s %8s %8s %8s %8s\n" "${SPECIES[$sp]}" "N/A" "N/A" "N/A" "N/A" "N/A"
    fi
done

echo ""
echo "Microhomology distribution (non-HR, all species combined):"
echo ""
printf "%6s %8s\n" "MH(bp)" "Count"
printf "%6s %8s\n" "------" "-----"

# Combine all TSVs and count microhomology for non-HR
for sp in ecoli_K12_MG1655 abaumannii_ATCC17978 scoelicolor_A3_2 mtuberculosis_H37Rv mpneumoniae_M129 paeruginosa_PAO1; do
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
