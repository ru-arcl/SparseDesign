#!/usr/bin/env bash
set -u

binary=${1:-./lineardesign-clean}
tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT
printf 'M\n' >"$tmpdir/protein"
tests=0

fail() {
    echo "CLI test failure: $*" >&2
    exit 1
}

expect_failure() {
    label=$1
    expected=$2
    shift 2
    tests=$((tests + 1))
    "$binary" "$@" <"$tmpdir/protein" >"$tmpdir/stdout" 2>"$tmpdir/stderr"
    status=$?
    if [ "$status" -eq 0 ]; then
        fail "$label unexpectedly succeeded"
    fi
    if ! grep -Fq -- "$expected" "$tmpdir/stderr"; then
        echo "stderr was:" >&2
        sed 's/^/  /' "$tmpdir/stderr" >&2
        fail "$label did not report '$expected'"
    fi
    if grep -Fq -- "terminate called" "$tmpdir/stderr"; then
        fail "$label terminated through an uncaught exception"
    fi
}

tests=$((tests + 1))
"$binary" --help >"$tmpdir/stdout" 2>"$tmpdir/stderr" || fail "--help failed"
grep -Fq -- "usage:" "$tmpdir/stderr" || fail "--help omitted usage"

tests=$((tests + 1))
"$binary" -l 0 -j 1 <"$tmpdir/protein" >"$tmpdir/stdout" 2>"$tmpdir/stderr" || \
    fail "valid options failed"
grep -Fq -- "mRNA sequence:" "$tmpdir/stdout" || fail "valid run omitted design output"

tests=$((tests + 1))
"$binary" --evaluate-mrna ATG <"$tmpdir/protein" >"$tmpdir/stdout" 2>"$tmpdir/stderr" || \
    fail "valid fixed-mRNA evaluation failed"
grep -Fq -- "mRNA sequence:  AUG" "$tmpdir/stdout" || \
    fail "fixed-mRNA evaluation did not normalize DNA to RNA"

expect_failure "missing -l value" "error: missing value for -l" -l
expect_failure "missing --codonusage value" "error: missing value for --codonusage" --codonusage
expect_failure "missing -j value" "error: missing value for -j" -j
expect_failure "missing motif value" "error: missing value for --forbid-motif" --forbid-motif
expect_failure "missing evaluated mRNA" "error: missing value for --evaluate-mrna" --evaluate-mrna
expect_failure "empty codon path" "requires a nonempty path" -c ""
expect_failure "empty motif" "requires a nonempty RNA motif" --forbid-motif ""
expect_failure "empty evaluated mRNA" "requires a nonempty RNA" --evaluate-mrna ""
expect_failure "invalid evaluated mRNA" "expects only A, C, G, U/T" --evaluate-mrna AXG
expect_failure "wrong evaluated mRNA length" "length does not equal three times" --evaluate-mrna AUGA
expect_failure "wrong evaluated translation" "does not encode the input protein" --evaluate-mrna UGG
expect_failure "incompatible evaluation option" "cannot be combined" \
    --evaluate-mrna AUG --sparse-stats
expect_failure "incompatible evaluation lambda" "cannot be combined" \
    --evaluate-mrna AUG -l 1

expect_failure "nonnumeric lambda" "--lambda expects a finite nonnegative number" -l nope
expect_failure "lambda suffix" "--lambda expects a finite nonnegative number" -l 1x
expect_failure "NaN lambda" "--lambda expects a finite nonnegative number" -l nan
expect_failure "infinite lambda" "--lambda expects a finite nonnegative number" -l inf
expect_failure "negative lambda" "--lambda expects a finite nonnegative number" -l -1
expect_failure "lambda that overflows internal scaling" \
    "--lambda expects a finite nonnegative number" -l 1e308

expect_failure "zero threads" "--threads expects a positive integer" -j 0
expect_failure "negative threads" "--threads expects a positive integer" -j -2
expect_failure "fractional threads" "--threads expects a positive integer" -j 1.5
expect_failure "thread suffix" "--threads expects a positive integer" -j 2x
expect_failure "oversized threads" "--threads expects a positive integer" -j 999999999999999999999

expect_failure "unknown long option" "error: unknown option: --unknown" --unknown
expect_failure "unexpected positional argument" "error: unknown option: protein.txt" protein.txt
expect_failure "missing codon table exception" "error: cannot open codon usage table:" \
    -c /definitely/not/a/codon-table.csv
expect_failure "invalid motif exception" "error: forbidden motifs must contain only" --forbid-motif AX
expect_failure "eliminating motif exception" "error: forbidden motifs eliminate every" --forbid-motif AUG

tests=$((tests + 1))
"$binary" --diagnostics <"$tmpdir/protein" >"$tmpdir/stdout" 2>"$tmpdir/stderr" || \
    fail "diagnostics failed"
grep -Fq -- "dp_cost=" "$tmpdir/stderr" || fail "diagnostics omitted DP objective"
grep -Fq -- "objective_delta=0" "$tmpdir/stderr" || fail "diagnostics reconstruction mismatch"

printf '>first\nM\n>second\nW\n' >"$tmpdir/protein"
expect_failure "multiple FASTA records" "expected one protein record"
printf 'M\n>second\nW\n' >"$tmpdir/protein"
expect_failure "mixed raw and FASTA" "expected one protein record"
printf '>first\n>second\nM\n' >"$tmpdir/protein"
expect_failure "empty first FASTA record" "expected one protein record"

echo "CLI parsing/error tests: $tests passed"
