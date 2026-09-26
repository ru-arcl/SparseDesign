# SparseDesign: exact codon-constrained RNA design

SparseDesign is an independently implemented C++ solver for a joint RNA folding-energy and
codon-usage objective over weighted codon automata. Candidate sparsification reduces multiloop
split work while preserving the specified dense recurrence; the complete solver remains
worst-case cubic in RNA length and quadratic in memory. See [DESIGN.md](DESIGN.md) for the
model, algorithm, numerical contract and implementation lineage.

[Benchmark reference results](BENCHMARKS.md) summarize the paper's primary AMD EPYC
measurements, with links to raw evidence, plots and reproduction instructions. The separate
commodity-PC example completes human Dp427c on a Core i9 workstation; its measurements
are not pooled with the EPYC results.

Original project code and documentation are **source-available under the
[Rutgers Non-commercial Research License (RU-NCRL)](LICENSE.md)**.
The license text is reproduced word for word from the supplied [Word document](LICENSE.docx).
See the license for the permitted Purpose, restrictions, derivative-work obligations and
Rutgers Technology Transfer contact. [NOTICE](NOTICE) carries the Rutgers copyright notice.
[Third-party notices](THIRD_PARTY_NOTICES.md) and [provenance](PROVENANCE.md) identify
separately licensed material. The restricted LinearDesign implementation and binaries are
not bundled.

## Requirements and build

The solver needs GNU Make and a C++11 compiler with OpenMP (tested with GCC on Linux).
The test harness additionally needs Bash and Python 3; analysis tools also use Python 3.
The commands below assume a Linux shell, including Linux under WSL on Windows. No ViennaRNA library,
Python package, network connection, or restricted binary is needed to build and run.
On Debian/Ubuntu the system build packages are `build-essential` and `python3`.
Clone the [SparseDesign repository](https://github.com/ru-arcl/SparseDesign), or copy this
`code/` directory outside the frozen publication bundle. Run the commands from the checkout's
root. The default codon table is relative to the
current directory. Use `-c /absolute/path/table.csv` when invoking the binary elsewhere.

```bash
git clone https://github.com/ru-arcl/SparseDesign.git
cd SparseDesign
make
printf ">demo\nMFLMVF\n" > protein.fasta
./sparsedesign -l 0.5 < protein.fasta
./sparsedesign -l 1.0 -c data/codon_usage_freq_table_human.csv < protein.fasta
./sparsedesign -l 1.0 -j 8 --sparse-stats < protein.fasta
./sparsedesign -l 1.0 --forbid-motif UUUUU --forbid-motif GGCGCC < protein.fasta
./sparsedesign --evaluate-mrna AUGUUCCUGAUGGUGUUC < protein.fasta
```

`sparsedesign` is the release executable. The build also provides the legacy name
`lineardesign-clean` for existing scripts; both run the same implementation. Historical
evaluation records retain the executable names and source identities used for those runs.

- `-l LAMBDA` — finite, nonnegative MFE–CAI tradeoff (0 = pure MFE; larger favors CAI). Default 0.
- `-c FILE` — codon-usage table (default human).
- `-j THREADS` — positive gap-wavefront worker count; default 1, capped by host topology and
  subject to OpenMP runtime limits. `-v` reports the effective count.
- `--forbid-motif RNA` — reject a motif anywhere in the CDS, including across codon boundaries;
  repeat the option for multiple motifs (`T` is accepted as `U`).
- `--sparse-stats` — emit publication/diagnostic candidate counts to standard error after an
  `O(num_nodes²)` read-only post-pass; reported storage covers candidate vectors, not total RSS.
- `--diagnostics` — emit full double-precision DP objective, reconstructed objective, their
  difference, codon penalty and CAI to stderr. Objective costs are in 0.01 kcal/mol units;
  the codon penalty `sum(-log(w))` and CAI are dimensionless. Every design
  checks its reconstructed weighted-DFA path and objective, even without this option.
- `--evaluate-mrna RNA` — validate that one supplied RNA/DNA sequence encodes the input protein,
  then report its Turner-d0 MFE structure and CAI without optimizing its codons. This mode
  rejects nonzero `--lambda`, `--forbid-motif`, and `--sparse-stats`.

Each invocation reads one protein in FASTA or raw form from stdin; multiple records are
rejected. It designs exactly the input residues: include `*` if a stop codon is required.
Standard output contains the input header (or `>designed`), labeled CDS and dot-bracket
structure lines, and a summary with MFE to two decimals and CAI to three decimals.
Use `--diagnostics` for full-precision objective comparisons; the displayed CAI is rounded.
Errors go to standard error and return exit status 1; successful runs return 0.

### Infer the supporting lambda range of an existing mRNA

An mRNA does not generally identify one lambda: it can be optimal over an interval, only at one
breakpoint, or at no nonnegative lambda. The inverse tool folds and scores the supplied synonymous
mRNA, adaptively searches the exact optimizer's supported MFE–CAI envelope, and reports a
numerical supporting interval:

```bash
python3 tools/frontier.py --infer-mrna target-mrna.fa -j 8 < protein.fa
python3 tools/frontier.py --infer-mrna target-mrna.fa --json < protein.fa
```

The tool locates `sparsedesign` and the bundled table relative to its own directory.
Use `--binary PATH` or `SPARSEDESIGN_BIN` to select another executable; the legacy
`LINEARDESIGN_CLEAN_BIN` environment variable remains supported as a fallback.

The upper endpoint search expands toward the minimum-codon-cost design, subject to explicit
iteration and lambda limits reported in JSON. This is a numerical supported-envelope calculation.
Unsupported nondominated designs are omitted, and ties share a representative. JSON distinguishes
coverage within numerical tolerances from exact completeness. An unbounded upper interval uses JSON
`"lambda_max": null` together with `"unbounded": true`. The calculation
can require roughly two exact design solves per supported frontier vertex and is therefore intended
for offline analysis. See the module docstring for tolerances and stopping conditions.

## Tests

```bash
make test          # optimized default-layout build and deterministic suite
make test-asan     # fresh AddressSanitizer binaries in /tmp
make test-ubsan    # fresh UndefinedBehaviorSanitizer binaries in /tmp
make lowmem-packed # explicit forward-only packed build for lower peak RSS
make test-layouts  # exact square/packed validation matrix
make test-lowmem-asan
make test-lowmem-ubsan
```

The exactness matrix includes exhaustive synonymous-CDS enumeration, dense/sparse weighted-lattice
equality, an exact-integer all-cell kernel, selected production multiloops, motif-product language
checks, dense-wavefront and scalar-right-dense equivalence at `j1`/`j2`/`j4`, all-cell
production Turner parity on 72 configurations, and CLI failure behavior. Acceptance is the
command's exit status, not a substring in its output.

The default production build uses the square memo arena, selected by the recorded Q6UXY8
layout screen at one and sixteen threads. `make lowmem-packed` builds the forward-only position-major arena
as `sparsedesign-lowmem-packed` (also available as `lineardesign-clean-lowmem-packed`). This
explicit memory-saving alternative retains the same exact recurrence and output contract,
with a runtime effect that depends on input and hardware. Both layouts are covered by `make test-layouts` and the
layout-specific sanitizer targets.

## Model and objective

The optimization minimizes `MFE_kcal + lambda * sum(-log(w(codon)))` over synonymous
CDSs and pseudoknot-free secondary structures. `w` is frequency divided by the maximum
synonymous frequency. The model is Turner 2004 at 37°C with dangles disabled (`d0`),
minimum hairpin loop size three, special tri/tetra/hexaloops enabled, and at most 30
unpaired nucleotides in a bulge/internal loop. Double arithmetic is used for weighted
objectives. Dense-equivalence guarantees apply under this model and in real arithmetic;
regression tests assess its floating-point implementation. MFE and CAI are computational
proxies and are not calibrated predictions of expression or therapeutic efficacy.

Custom codon tables use CSV rows `CODON,AMINO_ACID,FREQUENCY`, e.g. `UUC,F,0.54`.
Use one row per codon, uppercase one-letter amino-acid codes (`*` for stop), and finite
nonnegative frequencies; RNA or DNA spelling is accepted. The table defines the genetic
code used by the solver. Zero frequencies are floored to relative adaptiveness `1e-9`;
they do not exclude codons. Use motif constraints for exclusions. Include all desired
synonymous codons and all amino acids in the input. The bundled human and yeast tables
provide the data used by the manuscript's experiments.

## Controlled ablation and profiling APIs

`src/fold_turner.h` exposes three cost-only wavefront solvers:

- `turner_min_cost`: production candidate-sparse right-normal recurrence.
- `turner_min_cost_dense_right_normal`: the identical scalar recurrence and endpoint lists,
  retaining every finite direct interval. This isolates the retention predicate.
- `turner_min_cost_dense_wavefront`: the former pair-indexed left-normal branch scan.
  Its comparison measures the complete multiloop replacement.

The APIs accept lambda in **0.01 kcal/mol units**, so multiply the CLI lambda by 100.
`turner_min_cost_profiled` selects any kernel via `TurnerKernel` and records phase wall
times and split-list visits in `TurnerProfile`. Profiling uses separately compiled template
instantiations; ordinary solves perform no per-split instrumentation or clock reads.
Time ordinary kernels for performance comparisons and use the profiled path separately
for attribution. Inner-loop counters can substantially increase runtime; their phase times
do not quantify the ordinary solver's phase-time breakdown. `TurnerSparseStats` on the dense-right control reports all retained
finite direct intervals; those counts are not sparse candidates.

## Verify energy-parameter provenance

The source distribution contains generated headers. Regeneration is optional and requires
ViennaRNA's Python bindings with the pinned Turner-2004 parameter text (2.7.2 was used):

```bash
python3 -m venv .venv-vienna
.venv-vienna/bin/pip install ViennaRNA==2.7.2
python3 tools/generate_vienna_parameters.py --vienna-python .venv-vienna/bin/python --check
# Intentional regeneration, after checking provenance:
python3 tools/generate_vienna_parameters.py --vienna-python .venv-vienna/bin/python --write
```

`--check` fails on a changed upstream SHA-256 or any header difference except the source
version in the provenance comment. It never modifies the checked files. Parameter headers
are Makefile dependencies, so changing one rebuilds the executable.

## Verification scope

The default suite checks the production sparse and dense recurrences, reconstruction,
finite-state languages, small exhaustive optima and CLI failures using the default layout.
`make test-layouts` separately checks both storage layouts.
The optional independent validation campaign adds 130 prespecified cases, including published
regression fixtures, external ViennaRNA checks, full-precision objective reconstruction,
loop-cap boundaries and independent synonymous enumeration. Its [EPYC results and limitations](evaluation/epyc/validation/README.md)
identify which checks share solver energy functions and which use separate implementations.

[Historical timing evidence](evaluation/historical/README.md) preserves older binaries'
identifiers and separate protocols. These records do not assert that a newly compiled release
has the same runtime. Candidate counts describe retained inventory and cannot be interpreted
as whole-process speedups or memory savings. Long inputs can require multi-gigabyte working
sets; the CLI itself does not impose a time or memory limit.

## Layout
```
src/codon_table.*   genetic code, codon usage, w(c), CAI
src/dfa.*           weighted codon DFA (CAI -log w on edges)
src/energy.*        Turner 2004 d0 energy model (from open ViennaRNA/LinearFold formulas)
src/vienna/         generated Turner 2004 tables; see third-party notices
src/fold_simple.*   simplified Nussinov lattice DP (validation)
src/fold_turner.*   Turner-d0 lattice DP, backtrace and separate structure evaluation
src/main.cc         CLI
test/               unit/verification tests and design drivers
data/               codon-usage tables
tools/frontier.py   standalone numerical supported-envelope / inverse-lambda tool
tools/generate_vienna_parameters.py  optional pinned parameter verification/regeneration
```

## Publication evaluations

Start with [BENCHMARKS.md](BENCHMARKS.md) for the current reference results and their
measurement boundaries. The primary EPYC studies were completed on 25 September 2026.

The optional [external validation campaign](evaluation/epyc/validation/README.md) records independent
ViennaRNA 2.7.2 checks, external correctness fixtures, objective precision, length/composition
stress cases, and finite-state constraints. Its dependency is unnecessary for building or
running the solver. The [solver verification report](evaluation/solver/solver-report.md)
documents the dense-control APIs, production all-cell checks, layouts, sanitizer runs,
parameter-drift controls and packaged-tree checks. The [density study](evaluation/density/README.md)
reproduces the recorded 7,600-task analysis without another optimization campaign.
[Historical comparisons](evaluation/historical/README.md) keep the July protocols separate.
[Primary EPYC evidence](evaluation/epyc/README.md) contains the current three-arm kernel
comparison, separate profiling, complete-CLI scaling and Dp427c layout study, native
public-software comparison and validation rerun. All 934 non-pilot timing tasks succeed;
all 70 scaling/Dp427c outputs pass independent checks. The native comparison preserves
all feasibility and repeated outcomes, independently scored outputs and model-audit gates.
A separate 15-run EPYC follow-up is not pooled with the primary estimates.
The [commodity-PC dataset](evaluation/commodity-pc/README.md) retains only the ten i9
Dp427c runs used by the paper, with their inputs, measured sources and resource records.
The superseded i9 timing, native-comparison, validation and diagnostic campaigns are
omitted from this code release. [Release notes](RELEASE.md) describe the package scope
and verification.
For fresh timing measurements, set `export LC_ALL=C` before running collection commands;
the recorded GNU-time parsers expect English field names.

The source snapshot can be checked before building with `sha256sum -c SOURCE_SHA256SUMS`.
This manifest covers solver sources, bundled numerical tables, tools, tests and release
instructions; separate evaluation records identify the exact source used for each campaign.
`python3 tools/update_source_manifest.py --check` also checks for added or missing active source
files. After reviewing intentional changes, regenerate it with
`python3 tools/update_source_manifest.py`; this does not alter frozen evaluation manifests.
Before building, `python3 tools/verify_release.py` checks the complete source-and-evidence
package against `RELEASE_SHA256SUMS`, including the intended file set.

## Release maintenance

The standalone repository is [ru-arcl/SparseDesign](https://github.com/ru-arcl/SparseDesign).
This directory contains its complete buildable solver and publication evidence. For an
authorized release, keep `LICENSE.md`, `LICENSE.docx`, `NOTICE`, `THIRD_PARTY_NOTICES.md`,
`LICENSES/`, `CITATION.cff` and the evidence notices with the files they accompany. Run
`make test` and the checksum checks before releasing; keep citation metadata and
[benchmark references](BENCHMARKS.md) consistent with the released evidence. The `.gitignore`
excludes build outputs and Python environments. The [RU-NCRL](LICENSE.md) governs use and
redistribution; these maintenance instructions grant no additional permissions.
