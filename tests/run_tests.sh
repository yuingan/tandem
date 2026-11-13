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
for g in tandem_2x tandem_3x is_pair_30kb is_family_6x tandem_diverged tandem_short_150 \
         inverted_pair array_unit_113 array_partial_240; do
    out="$work/$g"
    tandem -r "$data/$g.fna" -o "$out" -t 1 > /dev/null 2>&1
    n_tandem=$(($(wc -l < "$out/tandem_duplications.tsv") - 1))
    n_disp=$(($(wc -l < "$out/dispersed_pairs.tsv") - 1))
    n_inv=$(($(wc -l < "$out/inverted_pairs.tsv") - 1))
    n_mult=$(($(wc -l < "$out/unit_multiples.tsv") - 1))
    n_short=$(($(wc -l < "$out/short_unit_arrays.tsv") - 1))
    printf "  %-18s %2d tandem  %2d dispersed  %2d inverted  %2d multiple  %2d short array\n" \
        "$g" "$n_tandem" "$n_disp" "$n_inv" "$n_mult" "$n_short"
done

echo
echo "Module 1: benchmark genome"
echo "--------------------------"
tandem -r "$data/benchmark.fna" -o "$work/benchmark" -t "${THREADS:-4}" > /dev/null 2>&1
python "$here/score_benchmark.py" "$data/benchmark_truth.json" "$work/benchmark"

status=0

echo
echo "Module 1: calls from self-alignments seen in real genomes"
echo "---------------------------------------------------------"
python "$here/check_module1_calls.py" || status=1

echo
echo "Module 1: copy sequences and report"
echo "-----------------------------------"
python "$here/check_report.py" "$work" "$data" || status=1

echo
echo "Microhomology at the junction"
echo "-----------------------------"
python "$here/check_microhomology.py" || status=1

echo
echo "Module 2: junction scan on simulated reads"
echo "-----------------------------------------"
python "$here/check_scan.py" "$work/scan" || status=1

echo
echo "Benchmark helpers"
echo "-----------------"
python "$here/breseq_junctions.py" check || status=1

exit $status
