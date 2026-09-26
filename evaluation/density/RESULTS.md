# Candidate-density findings

Verified all **7,600** recorded records across 11 merged products, including original CSV hashes, the canonical campaign-plan hash, exact task/panel membership, sequence identities, provenance fields and candidate/storage arithmetic. This work performs no new design optimization.

The two human panels contain 2,000 proteins each. Median candidate retention is 3.5282% at lambda 0 and 2.1544% at lambda 4; maxima are 43.7956% and 37.2263%. Median Z/n² is 0.01434125 and 0.00875663. All 2,000 paired proteins retain fewer candidates at lambda 4; the median ratio is 0.60917.

## Length sensitivity

| Minimum included aa | Proteins | Slope, lambda 0 | Slope, lambda 4 |
|---:|---:|---:|---:|
| 10 | 2,000 | 1.526702 | 1.487288 |
| 100 | 1,500 | 1.610286 | 1.586061 |
| 400 | 1,167 | 1.699236 | 1.677345 |
| 1000 | 747 | 1.793964 | 1.772580 |

The residual plot shows curvature. Continuous two-slope fits with a prespecified 400-aa knot give:

| Lambda | Lower slope | Upper slope | Reduction in log-space SSE |
|---:|---:|---:|---:|
| 0 | 1.423994 | 1.681581 | 24.06% |
| 4 | 1.367338 | 1.668169 | 27.52% |

This is descriptive model sensitivity, not an asymptotic complexity estimate. Full tables also include disjoint length-band fits and fixed knots at 100 and 1000 aa.

## Sequence-similarity sensitivity

Read-only MMseqs2 clustering yields 1,565 operational similarity groups: 1,328 singletons, maximum group size 20, and 24 groups crossing length strata. The recorded settings require at least 30% alignment identity and 80% coverage of both sequences, with additional prefilter/E-value settings. These are heuristic algorithm-defined clusters, not verified biological families.

Each interval below uses 10,000 paired replicates. Whole-cluster resampling preserves one multiplicity per cluster across lambdas and strata, then calibrates the original stratum totals. Equal-cluster fits give each cluster total weight one; their length-calibrated variant preserves the original stratum masses.

| Fit / bootstrap | Lambda 0: slope [95% interval] | Lambda 4: slope [95% interval] |
|---|---:|---:|
| Panel / within strata | 1.526702 [1.521009, 1.532318] | 1.487288 [1.481378, 1.493134] |
| Panel / whole clusters | 1.526702 [1.519873, 1.533457] | 1.487288 [1.479701, 1.494700] |
| Equal clusters | 1.518812 [1.511304, 1.526485] | 1.477424 [1.469284, 1.485653] |
| Equal clusters + length strata | 1.525177 [1.518205, 1.532105] | 1.485149 [1.477856, 1.492411] |

The length-range effect is larger than these cluster adjustments in this panel. Bootstrap intervals remain conditional sensitivity measures and cannot establish independent biological families or extrapolation validity.

## Codon tables and artifacts

The 10 codon tables use the same 400 proteins at lambda 4. Their point-fit slopes range from 1.376508 to 1.496665. Human and yeast frequencies are rounded; the other recorded tables retain full within-synonymous precision. Host fits report point estimates only.

The figures `density-sensitivity.pdf` and `host-comparison.pdf`, generated `claims.tex`, and LaTeX/CSV tables derive from the verified data. The raw bootstrap draws, cluster assignments, compressed original inputs and exact source snapshots are retained for reproduction. Candidate allocation describes candidate containers, not whole-process RSS; recorded instrumented wall time is not a controlled performance benchmark.
