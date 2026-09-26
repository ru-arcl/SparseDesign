# Solver and release engineering report

Prepared 2026-09-06/07 from `lineardesign-clean`. This file records the implementation and
verification used for the publication release; timed scientific results are recorded separately.

This is a historical engineering report. For the current package and benchmark
references, see [README](../../README.md) and [BENCHMARKS](../../BENCHMARKS.md).
Current scientific validation is under [epyc/validation](../epyc/validation/README.md),
and release terms are in [LICENSE.md](../../LICENSE.md). Source-manifest counts
and original workspace paths below describe the September 7 verification.

## Controlled ablation

The new `turner_min_cost_dense_right_normal` evaluates precisely the scalar right-normal
multiloop recurrence used by production: same pair-type minimum, endpoint-indexed vectors,
iteration order, memo layout, extension rules and wavefront. Its compile-time `Prune=false`
retention predicate stores every finite direct interval. The production `Prune=true` predicate
stores a direct interval only if its scalar cost is strictly below its partitionable/unpaired
alternative. Thus this comparison isolates the retention predicate and resulting list sizes.
The older dense-wavefront implementation is retained as a third arm, measuring the complete
multiloop replacement, including orientation and pair-type collapse.

Public signatures (internal lambda uses centikcal/mol, 100 times CLI lambda):

```cpp
double turner_min_cost_dense_right_normal(const DFA&, double lambda,
    int min_loop=3, int threads=1, TurnerSparseStats* = nullptr);

enum class TurnerKernel { Sparse, DenseRightNormal, DenseWavefront };
double turner_min_cost_profiled(const DFA&, double lambda, TurnerKernel,
    TurnerProfile&, int min_loop=3, int threads=1, TurnerSparseStats* = nullptr);
```

The profiled API uses separate template instantiations. Ordinary production solves execute
neither clocks nor per-split counters. `TurnerProfile` reports memo-construction wall time,
closed/gap phase wall time, multiloop phase wall time, external phase wall time, raw split-list
visits, and entries surviving endpoint-position guards. Phase measurements include scheduling
and barriers. Row counters have a unique writer at each wavefront diagonal; no atomics are used.
Profiled executions should diagnose work, while ordinary kernels establish performance.

## Correctness tests

`test/test_right_normal_cells.cc` compares every M1, M2 and closed-pair cost from the actual
production implementation against the original dense recurrence and the new dense-right
control. It covers 12 proteins (four targeted plus eight fixed-seed random), three lambda
values, two multiloop-unpaired penalties, and motif-product lattices: 72 configurations,
181,344 interval comparisons. Every configuration also checks that dense-right retains exactly
all feasible direct intervals and sparse retains a subset. Profiled objectives and counters
are checked at j1/j2/j4. Test-only cell accessors are absent from normal API builds.

The existing dense-wavefront tests now also verify the dense-right control at j1/j2/j4 on
true multiloops, nonzero multiloop penalties, and same-layer product-lattice states.

`make test` passed completely with GCC 14.3.0 and OpenMP. The suite includes exhaustive
synonymous-CDS enumeration, weighted-lattice equality, motif-language tests, true multiloops,
500 generated exact-integer abstract DAGs, the new production all-cell comparisons,
36 CLI error/output checks, and three inverse-lambda end-to-end cases.

Layout, sanitizer, second-compiler and packaged-tree checks all passed; final results are below.

## Release fixes

- CLI rejects multiple FASTA records and mixed raw/FASTA input rather than concatenating genes.
- `TurnerDesign` returns `realized_cost` and table-derived `codon_penalty`. Every design
  verifies that its reconstructed RNA follows the DFA and reaches its final state, then
  compares the DP optimum to independently evaluated structure energy plus realized DFA
  edge penalty. Acceptance tolerance is `1e-4 + 1e-10*max(abs(dp),abs(realized))` internal units.
- `--diagnostics` emits full 17-significant-digit DP/realized costs, objective delta, codon
  penalty and CAI to stderr. Standard output remains suitable for existing consumers.
- Portable `tools/frontier.py` contains the inverse-lambda/support-envelope tool with local
  defaults and only Python standard-library dependencies. The inverse-lambda test imports
  it and checks the API's JSON-safe `lambda_max=None` plus `unbounded=True` contract.
- Parameter generator enforces SHA-256
  `2a43345a495850cfd2e0a78c57c6e02085e6df3c53496fe3289dc294d21732ad`.
  Its real `--check` compares all four complete generated headers, ignoring only the
  upstream version string in the provenance comment. Both changed numeric data and a
  changed upstream source digest are verified to fail. No files are modified by checking.
- All parameter headers are Makefile dependencies.
- README documents build prerequisites, complete CLI examples, the 37°C/d0/30-unpaired
  internal-loop-cap model, numerical scope, custom codon tables, diagnostic/profiling
  semantics, optional parameter regeneration, inverse-lambda limits, and all test targets.
- DESIGN corrects the old ablation's attribution and documents the new controlled API.

## Verification records

- `solver-initial.log`: initial compile and dense/right/sparse objective tests.
- `solver-normal-tests.log`: complete optimized suite.
- `solver-parameter-verification.log`: valid pinned headers and both negative drift controls.
- `solver-test-layouts.log`: square/packed matrix.
- `solver-test-asan.log`, `solver-test-lowmem-asan.log`: AddressSanitizer layouts.
- `solver-test-ubsan.log`, `solver-test-lowmem-ubsan.log`: UndefinedBehaviorSanitizer layouts.

Compute-heavy checks acquire `compute.lock` so they do not overlap timed campaigns.

## Completed layout and relocation checks

`make test-layouts` passed with `LDCLEAN_VALIDATE_FORWARD_KEYS` under both square and
forward-packed arenas. Each layout includes the 72-case/181,344-interval production test,
full objective/language/multiloop tests, j1/j2/j4 equivalence and all 36 CLI checks.

An isolated temporary copy with no sibling web service passed the three inverse-lambda cases;
its frontend resolves both default executable and codon table inside the copy. A dry-run
Make dependency check verifies that modifying `src/vienna/intl22.h` schedules a complete
executable rebuild. See `solver-relocation-dependencies.log`.

Clang is not installed. GCC 15.2.0 is available as a second compiler version alongside the
primary GCC 14.3.0; its production-cell/CLI/inverse-lambda checks passed. This is a compiler
version check, not an independent compiler-family check.

## Standalone source package

`publication/code/` now contains source, local tools, tests, parameter/codon data, Makefile,
complete instructions, citation metadata and provenance notes. Binaries, caches and the
historical restricted-shared-library driver are excluded. The independent validation agent
owns `research/publication_validation.*` and `evaluation/validation/` inside that package.
No public license was inferred from another optimizer: the original clean tree did not
supply one, and selecting final release terms remains a rights-holder/release-maintainer task.

The standalone package includes `SOURCE_SHA256SUMS` covering 41 core source, test, tool,
numerical-data and instruction files. `sha256sum -c SOURCE_SHA256SUMS` passes. Scientific
validation inputs/results and their own source snapshots remain separately identified.

## Final verification matrix (2026-09-07)

| Check | Result | Record |
|---|---|---|
| Complete normal optimized suite, GCC 14.3 | PASS | `solver-normal-tests.log` |
| Square/packed with forward-key assertions | PASS | `solver-test-layouts.log` |
| AddressSanitizer, square | PASS | `solver-test-asan.log` |
| AddressSanitizer, packed | PASS | `solver-test-lowmem-asan.log` |
| UndefinedBehaviorSanitizer, square | PASS | `solver-test-ubsan.log` |
| UndefinedBehaviorSanitizer, packed | PASS | `solver-test-lowmem-ubsan.log` |
| GCC 15.2 actual-cell, CLI and inverse-lambda checks | PASS | `solver-gcc15.log` |
| Pinned parameter tables and negative drift controls | PASS | `solver-parameter-verification.log` |
| Relocated frontend and header rebuild dependency | PASS | `solver-relocation-dependencies.log` |
| Actual publication tree: fresh compile, CLI36, inverse3, parameter check | PASS | `solver-publication-smoke.log` |
| 41-file source manifest | PASS | `solver-source-checksums.log` |

AddressSanitizer used the existing harness setting `detect_leaks=0:halt_on_error=1`;
this checks address/lifetime errors but makes no leak-sanitizer claim. UBSan used
`halt_on_error=1:print_stacktrace=1`. Sanitizer targets rebuilt the complete deterministic
C++ and CLI suites in temporary directories, including all 181,344 production interval
comparisons under each layout. No sanitizer error was reported.

The source physics remained frozen throughout the contemporary scientific runs. The
validation agent separately reports all 130 external/precision/robustness/constraint
cases passing, plus independent enumeration and loop fixtures, with its exact snapshots
in `publication/code/evaluation/validation/`. No additional evaluation is needed to
establish the requested solver regression matrix.
