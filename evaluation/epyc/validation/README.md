# EPYC validation evidence

This is the manuscript's primary correctness and robustness validation dataset,
measured on `arrakis` (dual AMD EPYC 7313, 1 TiB RAM) on 2026-09-25. All 130
prespecified design cases, 28 loop fixtures and 30 exhaustive comparisons pass.
The 105-aa, 16-motif case at lambda 0 takes 37.6836 seconds and reaches
112,268 KiB (109.637 MiB) before dense checks. These are descriptive measurements
from one instrumented serial execution, not repeated performance estimates.

`designs.jsonl`, `exhaustive.json`, `exhaustive-sequences.jsonl`, `loops.json`,
`inputs.json`, `inputs/`, `build.json` and `validation.log` are copied unchanged
from `research/results/epyc-node1-2026-09/validation`. The collector's original,
less detailed summary is preserved as `collected-summary.json`. The complete
`summary.json` and `claims.tex` are rebuilt from these EPYC records with the
later frozen analysis script, without rerunning any solves. `raw-manifest.json`
records the hashes of every raw input and source file. `SHA256SUMS` covers the
whole dataset, excluding the checksum file itself.

The executed collector, C++ driver and solver sources are in `source/` and
match `build.json`; GCC 11.4.0 built the executable. ViennaRNA 2.7.2 independently
evaluated and refolded each returned design under the frozen Turner-d0 settings.
The later summary and plotting sources are separate in `source/analysis/`.
Comparison with the original i9 study confirmed that all 130 design records differ only in executable
hash and measured time/RSS: objectives, residuals, sequences, structures, motif
checks, candidate counts and independent folding results are identical. The
same comparison for exhaustive and loop records finds only time/RSS differences.
The superseded i9 validation dataset is omitted from this code release.

## Rebuild the manuscript evidence

From `publication/code`, Python 3.10+ and the standard library are sufficient
to verify source/raw hashes, recheck frozen case membership and translation,
and regenerate the complete summary and count macros:

```bash
python3 evaluation/epyc/validation/reproduce.py
```

To regenerate the figure and copy the summary, count macros and PDF into the
manuscript, use matplotlib (the recorded figure used version 3.10.9):

```bash
python3 evaluation/epyc/validation/reproduce.py --plot --paper ../paper
```

Without `--plot`, `--paper ../paper` imports the existing figure after checking
that its provenance matches the regenerated summary, raw designs and plotter.
All plotted candidate counts and automaton sizes agree exactly with the
original study; the figure is regenerated here to link it to the EPYC records.

## Check the numerical claims

The rebuilt summary reports 128 emitted designs and two expected infeasible
languages; all three dense controls agree on 114 cases. Its largest absolute
weighted-objective residual is `7.275957614183426e-11` centikcal/mol and the
largest external evaluation/refold difference is `2.4414062522737368e-5`
kcal/mol. Enumeration covers 1,931 feasible records (1,714 distinct RNAs),
including nine checks around three nonzero support-line breakpoints. The 28
loop fixtures include 22 special hairpins and six internal/bulge-loop boundaries;
the exactly constrained external recurrence accepts 30 unpaired bases and
rejects 31.

One of 15 full-precision/rounded-table pairs changes its RNA, with MFE increasing
by 4.60 kcal/mol and full-precision objective regret `0.008665995290444428`
kcal/mol. Seven of 96 strict-MFE within-protein comparisons reverse the AUP
ordering; neither d0/d2 MFE nor d0-MFE/d0-ensemble ordering reverses.
For 64/105-aa motif inputs, nodes grow from 217/362 to 598/966 and maximum
widths from 2/2 to 8/10. Lambda-0 constraint penalties are 0.50/3.40 kcal/mol,
and the largest weighted penalty is `8.716791318035611` kcal/mol. Every returned
RNA avoids its complete forbidden-motif inventory.

The release's production-cell, abstract-DAG, sanitizer and compiler checks in
the first paragraph of the manuscript's validation section are separate
software checks; this dataset supplies the external/stress panel and figure.

## Repeat the computation

To perform a new run, install `ViennaRNA==2.7.2`, copy `inputs.json` and `inputs/`
to a fresh output directory, and invoke the executed collector. For example,
from `publication/code`:

```bash
mkdir -p /tmp/sparsedesign-validation-repeat
cp evaluation/epyc/validation/inputs.json /tmp/sparsedesign-validation-repeat/
cp -R evaluation/epyc/validation/inputs /tmp/sparsedesign-validation-repeat/
python3 evaluation/epyc/validation/source/publication_validation.py \
  --source evaluation/epyc/validation/source/solver \
  --out /tmp/sparsedesign-validation-repeat --build \
  --lock /tmp/sparsedesign-validation-repeat/compute.lock
python3 evaluation/epyc/validation/source/analysis/publication_validation.py \
  --out /tmp/sparsedesign-validation-repeat --summarize-only
```

Elapsed time includes lattice construction, solving, traceback and candidate
statistics. Peak RSS is sampled before optional dense checks; `process_seconds`
includes those checks. The in-house evaluator and dense controls share energy
functions with the optimizer. Python codon-penalty checks and ViennaRNA provide
independent implementations. Returned-sequence refolding alone does not prove
global design optimality; the separate small exhaustive panel tests that claim.
