# SparseDesign EPYC source release — 26 September 2026

This package accompanies the technical report's primary AMD EPYC measurements.
[BENCHMARKS.md](BENCHMARKS.md) gives the reference numbers and measurement boundaries.
The active solver, numerical inputs and build behavior are unchanged by this update.
The supplied Rutgers license text and original Word document are unchanged.

This revision uses consistent terminology across the release guides and Markdown
benchmark reports, and includes the paper's direct Dp427c time/memory comparison.
The recorded analysis sources retain their measured identities and original wording;
regenerating a report may restore that wording while reproducing the same numerical results.
Dataset paths and schema keys remain stable for reproducibility.

## Included evidence

- [Primary EPYC studies](evaluation/epyc/README.md): 934 successful non-pilot timing
  runs, a separate 15-run follow-up, native-software comparisons, external validation
  and independent checks of all 70 complete-command timing outputs.
- [Commodity-PC example](evaluation/commodity-pc/README.md): the ten original i9
  Dp427c runs, with measured sources, inputs, raw outputs and a compact verifier.
- [Candidate density](evaluation/density/README.md),
  [historical comparisons](evaluation/historical/README.md) and
  [solver verification](evaluation/solver/solver-report.md).

The superseded full i9 timing, native-software, validation and diagnostic campaigns
and their development helpers are omitted. They remain in the private development
repository. The standalone repository starts with a fresh root commit for this release;
older standalone commit identifiers describe a superseded package.

The EPYC and density README paths were updated for this package layout, with their
containing manifests refreshed. Measured records, source snapshots, protocols and
numerical JSON/CSV results retain their original bytes. The compact commodity-PC dataset
identifies its subset explicitly and verifies each copied record against the original
study's identities. Historical reports describe the versions and licenses current
when their checks ran; the current terms are [LICENSE.md](LICENSE.md).

## Verify the release

Run these commands from the repository root before building:

```bash
python3 tools/update_source_manifest.py --check
python3 tools/verify_release.py
python3 evaluation/commodity-pc/verify.py
```

`SOURCE_SHA256SUMS` covers active code, tables, tests, tools and release instructions.
`RELEASE_SHA256SUMS` additionally covers every intended evidence file and its nested
manifest. The complete-package verifier checks both bytes and membership, excludes Git
metadata and Python caches/environments, and rejects native binaries and symlinks.
Use a separate working copy for builds, or run `make clean` before verifying again.
The supplied `.gitattributes` disables automatic line-ending conversion so checkout bytes
remain consistent with the manifests on every platform.
After reviewed changes, refresh the source manifest with `python3 tools/update_source_manifest.py`,
then the release manifest with `python3 tools/verify_release.py --write`.

The [README](README.md#tests) lists the solver tests. The EPYC, density and historical
datasets document analysis from frozen records; no new benchmark campaign is needed
to verify the reported results. Original absolute paths in collection logs describe
the measurement environment and are not required to analyze the release.

## Release checks

The unchanged solver's default `make test` suite passed on 25 September in a fresh,
isolated Linux directory with GCC 13.3.0. This includes 72 production-cell configurations (181,344 intervals), small
exhaustive and dense-control checks, 36 CLI cases and three inverse-lambda cases.
The packed executable also built and passed all 36 CLI and three inverse-lambda cases.
These are packaging checks, not new performance measurements. All nested evidence
manifests and active documentation links were checked. The trimmed evidence reproduces
the manuscript's numerical claims and five figures without changing their values.
