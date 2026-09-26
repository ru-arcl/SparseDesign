#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/ldclean-layouts.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT

cxx=${CXX:-g++}
common=(-std=c++11 -O2 -Wall -fopenmp)
solver_sources=(
    "$root/src/fold_turner.cc"
    "$root/src/energy.cc"
    "$root/src/dfa.cc"
    "$root/src/codon_table.cc"
)
tests=(test_sparse_exact test_motif_dfa test_turner_multiloop test_dense_wavefront)

run_layout() {
    local layout=$1
    shift
    local flags=("$@")
    for test_name in "${tests[@]}"; do
        "$cxx" "${common[@]}" "${flags[@]}" \
            "$root/test/${test_name}.cc" "${solver_sources[@]}" \
            -o "$tmpdir/${layout}-${test_name}"
        (cd "$root" && "$tmpdir/${layout}-${test_name}")
    done

    "$cxx" "${common[@]}" "${flags[@]}" \
        "$root/test/test_right_normal_cells.cc" "$root/src/energy.cc" \
        "$root/src/dfa.cc" "$root/src/codon_table.cc" -o "$tmpdir/${layout}-right-cells"
    (cd "$root" && "$tmpdir/${layout}-right-cells")

    "$cxx" "${common[@]}" "${flags[@]}" \
        "$root/src/main.cc" "${solver_sources[@]}" \
        -o "$tmpdir/${layout}-lineardesign-clean"
    (cd "$root" && bash test/test_cli.sh "$tmpdir/${layout}-lineardesign-clean")
}

run_layout square -DLDCLEAN_SQUARE_MEMOS -DLDCLEAN_VALIDATE_FORWARD_KEYS
run_layout packed -DLDCLEAN_VALIDATE_FORWARD_KEYS

echo "layout matrix: PASS (square, packed)"
