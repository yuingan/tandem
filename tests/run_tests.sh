#!/usr/bin/env bash
# Build test genomes and reads, run Modules 1 and 2 over them, check answers.
#
#   bash tests/run_tests.sh
#
# Needs tandem on PATH and MUMmer installed. Everything generated goes under
# $TANDEM_TEST_DIR (default /tmp/tandem_tests), nothing inside the repo.

set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
work=${TANDEM_TEST_DIR:-/tmp/tandem_tests}
data="$work/data"

mkdir -p "$data"
python "$here/make_test_genomes.py" "$data" > /dev/null
python "$here/make_benchmark_genome.py" "$data" > /dev/null

echo "Module 1: single-feature genomes"
echo "--------------------------------"
for g in tandem_2x tandem_3x is_pair_30kb is_family_6x tandem_diverged tandem_short_150; do
    out="$work/$g"
    tandem -r "$data/$g.fna" -o "$out" -t 1 > /dev/null 2>&1
    n_tandem=$(($(wc -l < "$out/tandem_duplications.tsv") - 1))
    n_disp=0
    [ -f "$out/dispersed_pairs.tsv" ] && n_disp=$(($(wc -l < "$out/dispersed_pairs.tsv") - 1))
    printf "  %-18s %2d tandem  %2d dispersed\n" "$g" "$n_tandem" "$n_disp"
done

echo
echo "Module 1: benchmark genome"
echo "--------------------------"
tandem -r "$data/benchmark.fna" -o "$work/benchmark" -t "${THREADS:-4}" > /dev/null 2>&1
python "$here/score_benchmark.py" "$data/benchmark_truth.json" "$work/benchmark"

echo
echo "Module 2: junction scan on simulated reads"
echo "-----------------------------------------"
python "$here/check_scan.py" "$work/scan"
