# Primary EPYC evidence — 25 September 2026

The revised paper uses the dual-16-core AMD EPYC 7313 server (1 TiB RAM,
GCC 11.4.0) for its primary timing, native-baseline and validation measurements.
Benchmarks ran on a shared server, introducing variability between runs, although
overall load was generally not exceptionally high when runs were attempted.
The separate [commodity-PC dataset](../commodity-pc/README.md) retains the ten
Intel Core i9-14900KF Dp427c runs used by the revised manuscript. Other superseded
i9 campaigns are omitted from this code release; platforms are never pooled.

| Study | Completed evidence |
| --- | --- |
| `timings/` | 720 ablation, 144 profile, 60 scaling and 10 Dp427c runs, all successful |
| `followup-Q8VIM6/` | 15 successful runs in five separate three-kernel blocks |
| `public-baselines/` | 96 feasibility outcomes and 210 successful fresh repetitions; parameter audit and independent output checks |
| `validation/` | 130 design cases, 28 loop fixtures and 30 exhaustive comparisons pass |
| `cli-validation/` | All 70 scaling/Dp427c outputs independently checked; five distinct RNA/structure pairs |
| `load-audit.json` | Per-run mean CPU utilization and exact raw-record hashes |

Five excluded pilots were intentionally not collected. Nineteen pilot records
are retained, including one timeout, and are excluded from performance estimates.
Thus the timing analyzer's overall `complete` flag is false even though all 934
non-pilot tasks are complete. The ordinary wall cap is 1,800 s, the Dp427c cap
3,600 s, and the native comparison cap 240 s. Timing/native measurements retain
the 32-GiB address-space cap. This cap is distinct from measured peak RSS and
from physical host memory. Frozen build records contain exact compiler flags,
source identities and machine descriptions.

`load-audit.json` retains the available per-run CPU snapshots and their coverage
limits as diagnostic evidence; the paper makes no numerical utilization guarantee.

## Reproduce without solving

From this directory (Python 3; plotting also requires matplotlib):

```bash
sha256sum -c SHA256SUMS
python3 -B timings/analyze_node1.py --input timings --output /tmp/epyc-timing-analysis --plot
python3 -B analyze_followup.py
python3 -B audit_cpu_utilization.py --input . --output /tmp/epyc-load-audit.json
python3 -B validate_cli.py --input timings --output /tmp/epyc-cli-check --lock /tmp/epyc-cli-check.lock --allow-partial
```

The CLI check needs ViennaRNA 2.7.2 and only evaluates returned structures; it
performs no design or refolding. `--allow-partial` accommodates the five missing
pilots; the manuscript importer separately requires all 70 CLI checks to pass.
The CLI summary's `original_study_complete: false` refers to those pilots.
See the native-comparison and validation subdirectory instructions for their
independent analysis dependencies and commands. Native executables, libraries
and license-restricted upstream sources are not distributed.

`timings/source/`, the build/protocol records and raw outputs preserve the
measured source identities. `collection/` records the execution drivers;
`timings/rejected-runs/` preserves replaced observations and its dataset index.
Only `timings/timings/runs/` contributes to primary estimates. Original commands
and absolute paths are provenance, not a requirement to reproduce analysis.
The EPYC Q8VIM6 follow-up uses the original registered input/order but is a
separate EPYC experiment; it makes no quiet-host claim and is not pooled.

From `publication/paper/`, `python3 -B tools/import_timing_evidence.py` now
imports the EPYC timing and native studies by default, preserving source hashes
in the paper's data directory. The older import path requires `--legacy-i9`.
The validation reproduction helper imports its own EPYC summary/claims/figure.
Re-importing data changes no solver or statistical estimator.
