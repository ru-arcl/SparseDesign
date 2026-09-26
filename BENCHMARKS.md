# SparseDesign benchmark reference results

These results accompany the revised SparseDesign technical report. The primary performance
and external-validation studies were completed on 25 September 2026 on a dual-16-core
AMD EPYC 7313 server with 1 TiB RAM, using GCC 11.4.0. Benchmarks ran on a shared server,
introducing variability between runs, although overall load was generally not exceptionally
high when runs were attempted. Results identify the measured frozen sources; recompiling
the release on another system does not guarantee the same timings.

The [EPYC evidence index](evaluation/epyc/README.md) links the protocols, source identities,
raw outputs, integrity checks and analysis commands. All **934 non-pilot timing runs**
succeeded: 720 ablation, 144 profile, 60 scaling and 10 Dp427c runs. Five excluded pilots
are intentionally absent, so the analyzer's overall `complete` flag remains false.
The separate 15-run Q8VIM6 follow-up is not pooled with these estimates.

## Controlled kernel comparison

The panel contains 12 selected natural proteins, each measured in five paired repetitions
at each lambda/thread setting. Ratios below are geometric means of the condition-wise
median paired cost-call ratios. Values above one favor the sparse kernel.

| Lambda | Threads | Scalar dense / sparse | Typed dense / sparse |
| --- | ---: | ---: | ---: |
| 0 | 1 | 1.055 | 1.321 |
| 0 | 16 | 1.091 | 1.399 |
| 4 | 1 | 1.067 | 1.326 |
| 4 | 16 | 1.114 | 1.424 |

The scalar dense control differs from sparse only in its retention predicate, isolating
pruning and the resulting list sizes. The typed dense comparison measures the complete
multiloop replacement, including orientation, scalar collapse and branch representation.
Cost-call time includes solver construction, forward evaluation and destruction; it
excludes input/DFA construction and traceback. Aggregate ratios do not imply a speedup
for every protein. These are same-source controls, not upstream LinearDesign comparisons.

See the [aggregate CSV with conditional 95% intervals](evaluation/epyc/timings/analysis/ablation-aggregate.csv),
[per-protein plot](evaluation/epyc/timings/analysis/controlled-ablation.png) and
[measurement definitions](evaluation/epyc/timings/analysis/METHODS.md).
The bootstrap intervals describe repeat variability within this fixed panel, not
population uncertainty or systematic hardware drift. Profiling uses separate instrumented
runs; its counters and phase times are not whole-command speedups.

## Complete-command scaling and memory

For Q54TT4 and Q61879 at lambda 0 and 4, ratios of three-run median command times give
**10.00–12.77-fold speedup from one to 16 threads**. These measurements include startup,
DFA construction, design, traceback and output. The full thread ladder is in
[scaling.csv](evaluation/epyc/timings/analysis/scaling.csv).

Human Dp427c (RefSeq NP_000100.3; 3,677 amino acids, 11,031 nucleotides) was measured
at lambda 0 with 16 threads. Each layout completed five runs:

| EPYC layout | Median command time (s) | Median peak RSS (GiB) |
| --- | ---: | ---: |
| Square (default) | 237.39 | 20.09 |
| Packed (`make lowmem-packed`) | 236.54 | 14.43 |

The median paired RSS reduction is **28.20%**; the median paired packed/square time ratio
is **0.996**. This layout comparison does not isolate candidate pruning. See the
[Dp427c summaries](evaluation/epyc/timings/analysis/workstation.csv) and
[full timing analysis](evaluation/epyc/timings/analysis/summary.json).
All 70 scaling/Dp427c outputs pass the recorded
[independent CLI checks](evaluation/epyc/cli-validation/REPORT.md).

On the same EPYC 7313 server, the single-thread local dense LinearDesign fork takes
4,912 seconds and 402.10 GiB peak RSS for this input at lambda 0. Relative to that
measurement, 16-thread packed SparseDesign delivers a **20.8-fold wall-clock speedup**
and a **27.9-fold reduction in peak memory use**. The ratios divide the single-thread
observation by the respective five-run packed medians; they include implementation
and thread-count differences. The comparator's `ld_j1` record is retained in the
[large-instance results](evaluation/historical/summary.json).

Timing runs used a 32-GiB address-space cap, with a 1,800-second ordinary wall cap and
a 3,600-second Dp427c wall cap. An address-space cap is distinct from measured peak RSS
and does not simulate a physically 32-GiB machine. The standalone CLI does not impose
these resource limits itself.

## Native software and validation

The [native software comparison](evaluation/epyc/public-baselines/README.md) measures
pinned LinearDesign, LinearCDSfold, DERNA and SparseDesign with one thread, a 240-second
wall cap and a 32-GiB address-space cap. It retains 96 feasibility outcomes: 42 successes,
8 timeouts and 46 prescribed skips. All 210 fresh repetitions succeeded. Skipped conditions
remain unattempted and do not establish maximum feasible input lengths.

Only four lambda-0 DERNA/SparseDesign conditions pass the thermodynamic-parameter,
admissibility, energy-agreement and output-validation gates. The DERNA/SparseDesign
median command-time ratios are **5.08, 4.71, 5.51 and 6.04** for **82, 95, 255 and 310
amino acids**, respectively. See the
[matched-ratio data](evaluation/epyc/public-baselines/analysis/matched-lambda0-ratios.json)
and [native time/memory plot](evaluation/epyc/public-baselines/analysis/native-cli-time-memory.png).
Native models otherwise differ; the paper excludes strict LinearDesign/LinearCDSfold
and lambda-4 speed ratios. Validating returned outputs alone does not establish
common-model design optimality.

The [external validation dataset](evaluation/epyc/validation/README.md) records
**130/130 design cases, 28/28 loop fixtures and 30/30 exhaustive comparisons passing**.
The 105-amino-acid case with 16 forbidden motifs takes **37.7 s and 109.6 MiB** in one
instrumented serial run, including construction, solving, traceback and statistics;
this is a descriptive observation, not a repeated timing estimate.

## Separate commodity-PC demonstration

On a Core i9-14900KF workstation with approximately 64 GiB RAM and GCC 14.3.0, the
same Dp427c input at lambda 0 with 16 threads completes five runs per layout:

| Core i9 layout | Median command time (s) | Median peak RSS (GiB) |
| --- | ---: | ---: |
| Square (default) | 128.26 | 20.10 |
| Packed | 126.42 | 14.43 |

These [recorded workstation results](evaluation/commodity-pc/analysis/workstation.csv) demonstrate
commodity-PC feasibility. They use a 32-GiB address-space cap and a 1,800-second wall cap;
the 16 threads span performance and efficiency cores. They remain separate from the
primary EPYC estimates. Older July comparisons also remain separate in the
[historical evidence dataset](evaluation/historical/README.md).

## Candidate counts and reproduction

The earlier [7,600-task candidate-density study](evaluation/density/README.md) is unchanged.
Across the 2,000-protein natural panel, median retained fractions are 3.53% at lambda 0
and 2.15% at lambda 4, approximately 28-fold and 46-fold inventory reductions.
Synthetic stress families can benefit little. Candidate inventory reductions are neither
whole-command speedups nor whole-process memory reductions.

Follow the [EPYC analysis instructions](evaluation/epyc/README.md#reproduce-without-solving)
to verify hashes and regenerate summaries and plots from frozen records without rerunning
designs. The [native analysis](evaluation/epyc/public-baselines/README.md#verify-and-reproduce-the-analysis)
and [validation reproduction](evaluation/epyc/validation/README.md#rebuild-the-manuscript-evidence)
describe their additional dependencies. Frozen protocols and source snapshots identify
the exact inputs and implementations; ordinary release builds and tests are documented
in [README.md](README.md).
