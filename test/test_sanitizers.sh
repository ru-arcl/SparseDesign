#!/usr/bin/env bash
set -euo pipefail

mode=${1:-}
case "$mode" in
    asan|ubsan) ;;
    *)
        echo "usage: $0 {asan|ubsan} [extra compiler flags ...]" >&2
        exit 2
        ;;
esac
shift

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/ldclean-${mode}.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT

cxx=${CXX:-g++}
flags=(-std=c++11 -O1 -g -Wall -fopenmp -fno-omit-frame-pointer)
if [[ $mode == asan ]]; then
    flags+=(-fsanitize=address)
    # LeakSanitizer is unreliable under the ptrace-based execution environment used for the
    # historical checks. Address, use-after-free, bounds, and lifetime checks remain enabled.
    export ASAN_OPTIONS=${ASAN_OPTIONS:-detect_leaks=0:halt_on_error=1}
else
    flags+=(-fsanitize=undefined)
    export UBSAN_OPTIONS=${UBSAN_OPTIONS:-halt_on_error=1:print_stacktrace=1}
fi
flags+=("$@")

build_and_run() {
    local name=$1
    shift
    "$cxx" "${flags[@]}" "$@" -o "$tmpdir/$name"
    (cd "$root" && "$tmpdir/$name")
}

build_and_run test_codon \
    "$root/test/test_codon.cc" "$root/src/codon_table.cc"
build_and_run test_dfa \
    "$root/test/test_dfa.cc" "$root/src/dfa.cc" "$root/src/codon_table.cc"
build_and_run test_fold_simple \
    "$root/test/test_fold_simple.cc" "$root/src/fold_simple.cc" \
    "$root/src/dfa.cc" "$root/src/codon_table.cc"
build_and_run test_energy_threadsafe \
    "$root/test/test_energy_threadsafe.cc" "$root/src/energy.cc" -pthread
build_and_run test_sparse_exact \
    "$root/test/test_sparse_exact.cc" "$root/src/fold_turner.cc" \
    "$root/src/energy.cc" "$root/src/dfa.cc" "$root/src/codon_table.cc"
build_and_run test_motif_dfa \
    "$root/test/test_motif_dfa.cc" "$root/src/fold_turner.cc" \
    "$root/src/energy.cc" "$root/src/dfa.cc" "$root/src/codon_table.cc"
build_and_run test_turner_multiloop \
    "$root/test/test_turner_multiloop.cc" "$root/src/fold_turner.cc" \
    "$root/src/energy.cc" "$root/src/dfa.cc" "$root/src/codon_table.cc"
build_and_run test_dense_wavefront \
    "$root/test/test_dense_wavefront.cc" "$root/src/fold_turner.cc" \
    "$root/src/energy.cc" "$root/src/dfa.cc" "$root/src/codon_table.cc"
build_and_run test_sparse_kernel "$root/test/test_sparse_kernel.cc"
build_and_run test_right_normal_cells "$root/test/test_right_normal_cells.cc" \
    "$root/src/energy.cc" "$root/src/dfa.cc" "$root/src/codon_table.cc"

"$cxx" "${flags[@]}" \
    "$root/src/main.cc" "$root/src/fold_turner.cc" "$root/src/energy.cc" \
    "$root/src/dfa.cc" "$root/src/codon_table.cc" -o "$tmpdir/lineardesign-clean"
(cd "$root" && bash test/test_cli.sh "$tmpdir/lineardesign-clean")

echo "$mode suite: PASS"
