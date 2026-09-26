# Recorded candidate-count analysis

This artifact analyzes existing results from the frozen July 2026 campaign.
It does not call a design optimizer, produce new coding sequences or alter the
recorded measurements. Generated findings are in `RESULTS.md`; `summary.json`
contains full-precision statistics. `claims.tex` and `tables/*.tex` supply
paper-ready values without manually entered numerical results.

## Reproduce from the standalone code repository

No external Python libraries, MMseqs2 installation or development checkout
are needed to verify all inputs and regenerate summaries from the frozen draws:

```bash
python3 evaluation/density/source/publication_density.py --out evaluation/density
python3 evaluation/density/source/test_publication_density.py --out evaluation/density
```

The first command verifies original decompressed CSV hashes, the original
canonical campaign-plan hash, exact planned membership, source/table identities,
ratio/allocation arithmetic, and the clustering input/partition provenance.
It regenerates every JSON, CSV and LaTeX summary. The second command checks the
earlier independent campaign audit, weighted regression against Python's
independently implemented regression, paired bootstrap scale invariance, changed
CSV rejection and invalid cluster-membership rejection.

To regenerate the actual Monte Carlo draws with the standard library:

```bash
python3 evaluation/density/source/publication_density.py --out evaluation/density \
  --rebootstrap --replicates 10000 --lock /tmp/sparsedesign-evaluation.lock
```

Use the same lock path as any concurrently running evaluation. The default
summary command reuses and verifies the frozen bootstrap cache. The cache is
bound to the recorded products, cluster assignments, seeds, replicate count
and bootstrap schema; input drift causes an error.

For the optional publication figures:

```bash
python3 -m venv /tmp/sparsedesign-figures
/tmp/sparsedesign-figures/bin/pip install matplotlib==3.10.6
/tmp/sparsedesign-figures/bin/python evaluation/density/source/publication_density.py \
  --out evaluation/density --plot
```

The two vector PDFs have inspected PNG previews. Their provenance records the
plotting package version, generator hash and summary hash. Scatter points in the
residual panel are rasterized inside the PDF; text and linework remain vector.

## Statistical definitions

The estimand is the observed, deliberately length-balanced protein panel. These
are **not** population-weighted descriptions of Swiss-Prot. The two human-table
objectives have 2,000 proteins each; ten codon tables share a 400-protein subset
at lambda 4, with the human subset reused from its larger panel. Quantiles use
linear interpolation between adjacent ordered observations.

`Z` counts retained branch candidates, `direct` counts feasible scalar direct
intervals, and `n` is RNA length (three times protein length). `Z/direct`
describes inventory contraction; it does not count executed split evaluations.
Candidate allocation obeys the recorded implementation's `24*N + 16*capacity`
formula and is not total RSS.

Global and cutoff fits are unweighted OLS of `ln(Z)` on `ln(n)`. Disjoint fits
use 10–99, 100–399, 400–999 and 1000–3999-aa bands. Separate continuous hinge
fits use each prespecified knot at 100, 400 or 1000 aa:

```
ln(Z) = intercept + lower_slope * ln(n)
        + slope_change * max(0, ln(n) - ln(3*knot_aa))
```

The upper slope is `lower_slope + slope_change`. SSE reductions are descriptive
comparisons with the single-slope fit, without a significance or complexity
claim. Residuals are retained per observation and summarized per length stratum.

Each bootstrap scheme uses 10,000 paired replicates and no finite-population
correction. The within-stratum sequence bootstrap uses seed 20260716; whole-
cluster sampling uses a separate RNG with seed 20260717. Both lambda endpoints
use identical sampled sequence/cluster multiplicities in each replicate.

Four results are reported separately:

1. The original unweighted fit with within-stratum sequence resampling.
2. The same fit with whole similarity clusters resampled uniformly, preserving
   their members across all strata, followed by calibration of each resampled
   stratum back to its original sequence count.
3. A different descriptive estimand giving each cluster total weight one
   (`1/cluster_size` per sequence), with whole-cluster bootstrap intervals.
4. Equal-cluster weights additionally calibrated within length strata to retain
   the original stratum masses, with the analogous whole-cluster bootstrap.

This distinguishes changes in dependence assumptions from changes in weighting
and length composition. In particular, equal-cluster weighting alone need not
preserve the original balanced length distribution. No resample may silently
omit an entire stratum; the implementation fails if that occurs.

The MMseqs2 grouping is a sensitivity proxy. It neither identifies all homologs
nor establishes independent biological families. Consequently, even the cluster
intervals are conditional diagnostics, not general biological confidence
statements. Codon-table comparisons retain the original point-fit-only convention
on their shared subset.

## Clustering provenance and optional rerun

`clustering/` retains the exact accession-sorted input FASTA (compressed), an
index with sequence hashes and length strata, raw representative/member
assignments, full command log, `/usr/bin/time -v` output, source-input hashes,
binary hash and recorded MMseqs2 version. The clustering binary is not distributed.

The recorded MMseqs2 version string is
`eec9c354be4276d2373996af2e50808b1390d527`. The one-step greedy length clustering
uses actual alignment identity >=0.30 (`--alignment-mode 3 --seq-id-mode 0`),
coverage >=0.80 of both sequences (`--cov-mode 0`), E-value <=0.001, sensitivity
7.5, at most 3,000 prefilter hits per query, four CPU threads and a 4-GiB split
memory limit. Default low-complexity masking and composition-bias correction
remain enabled. The group is an algorithm output at these settings; membership
does not imply that every pair in the group passes a direct pairwise test.
See the primary [MMseqs2 repository and documentation](https://github.com/soedinglab/MMseqs2)
and [parameter definitions](https://github.com/soedinglab/MMseqs2/blob/master/src/commons/Parameters.cpp).

To rerun the read-only clustering in a fresh directory with an installed
MMseqs2 binary:

```bash
mkdir -p /tmp/sparsedesign-recluster
cp evaluation/density/clustering/panel.fasta.gz /tmp/sparsedesign-recluster/
cp evaluation/density/clustering/sequence-index.csv /tmp/sparsedesign-recluster/
cp evaluation/density/clustering/input-provenance.json /tmp/sparsedesign-recluster/
python3 evaluation/density/source/publication_density_cluster.py \
  --out /tmp/sparsedesign-recluster --mmseqs /path/to/mmseqs \
  --lock /tmp/sparsedesign-evaluation.lock
```

The helper refuses to overwrite existing assignments and records the actual
binary and command. A different software version or clustering output is a new
sensitivity run and requires newly generated bootstrap draws. No optimizer is
part of this procedure.

## Files and attribution

`inputs/` is a compact, lossless bundle of the eleven original merged CSVs, both
panel manifests, the campaign plan and merged manifest. Gzip compression uses
a zero timestamp; decompressing reproduces the recorded bytes and their original
SHA-256 hashes. Legacy absolute paths remain in the original provenance, but
the analyzer resolves only validated local product names. Original panel
manifests preserve the UniProt/Swiss-Prot 2026_02 provenance and source metadata.
`clustering/panel.fasta.gz` is derived from the frozen original panel ZIP, verified
against every individual sequence and FASTA hash; it is used only for read-only
similarity analysis.

`bootstrap-draws.csv.gz` retains every draw; `bootstrap-provenance.json` records
its hash and executed generator version. `source/bootstrap-publication_density.py`
preserves that exact initial generator. The current source snapshot adds summary
and table formatting while using the same frozen draws. `verification/` contains
the earlier independent audit and the new numerical/integrity check log.
`SHA256SUMS` covers the final artifact. Neither external tool code nor the
restricted historical baseline binaries are included.

The original research tree can regenerate the compact bundle with
`research/publication_density.py --prepare --out PATH`, after copying the frozen
`clustering/` directory to that output. Input export reads the original datasets;
it never regenerates their candidate-count measurements.
