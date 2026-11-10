#!/usr/bin/env bash
# Build the test genomes, run Module 1 over them, and check the answers.
#
#   bash tests/run_tests.sh
#
# Needs tandem on PATH and MUMmer installed.

set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
data="$here/data"
work=${TANDEM_TEST_DIR:-/tmp/tandem_tests}

mkdir -p "$work"
python "$here/make_test_genomes.py" "$data" > /dev/null
python "$here/make_benchmark_genome.py" "$data" > /dev/null

echo "Single-feature genomes"
echo "----------------------"
for g in tandem_2x tandem_3x is_pair_30kb is_family_6x tandem_diverged tandem_short_150; do
    out="$work/$g"
    tandem -r "$data/$g.fna" -o "$out" -t 1 > /dev/null 2>&1
    n_tandem=$(($(wc -l < "$out/tandem_duplications.tsv") - 1))
    n_disp=0
    [ -f "$out/dispersed_pairs.tsv" ] && n_disp=$(($(wc -l < "$out/dispersed_pairs.tsv") - 1))
    printf "  %-18s %2d tandem  %2d dispersed\n" "$g" "$n_tandem" "$n_disp"
done

echo
echo "Benchmark genome"
echo "----------------"
tandem -r "$data/benchmark.fna" -o "$work/benchmark" -t "${THREADS:-4}" > /dev/null 2>&1
python "$here/score_benchmark.py" "$data/benchmark_truth.json" "$work/benchmark"
