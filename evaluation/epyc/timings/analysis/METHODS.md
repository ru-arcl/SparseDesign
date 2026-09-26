**cost seconds.** Harness steady-clock interval around cost call: solver construction, fill and destruction; excludes table/DFA construction, effective-team probe and traceback.

**cli wall seconds.** GNU-time complete child process: startup, DFA, design, traceback, output and destruction. Python Popen/wait wall is retained separately.

**profile seconds.** Separate instrumented runs. Construction is solver-constructor elapsed time. Closed, multiloop and external sum wavefront phase elapsed timers including barriers; not CPU-seconds. Dense typed closed phase also builds eager branch rows. Phase sum excludes teardown/counter setup/collection and other call overhead.

**split visits.** List elements entered in multiloop scans, including entries rejected by endpoint/start checks. Typed scans can include a terminating future-end entry. Not candidate count Z.

**split eligible.** Scalar/sparse: visited entries with start layer beyond interval start. Typed: excludes direct/same-layer/future-end cases. Pruning comparison uses scalar/sparse in the same right-normal loop.

**visit ratios.** Dense-scalar counter divided by sparse counter for a matching profile case; null when denominator is zero. Not wall speedup or candidate retention fraction.

**peak rss gib.** GNU maximum resident set size KiB / 2^20; whole-process memory, distinct from the 32GiB virtual address-space limit.

**swapping.** Host-wide /proc/vmstat pswpin/pswpout differences. Positive values cannot be attributed to this child; zeros show no host swapping during that window.

**partial outcomes.** Missing, timeout, failed and invalid-output outcomes remain in every planned grid. Timeouts are censored elapsed observations, not completed solve times. Successful-subset memory/timing summaries retain explicit counts.

**provenance.** Source/input/output hashes and recorded commands are checked; present binaries are rehashed and absent binaries may be omitted for redistribution. Actual inherited environment, affinity success and exclusive hardware use are not independently authenticated by these records.

Bootstrap specification:

```json
{
  "draws": 10000,
  "seed": 20260908,
  "estimand": "geometric mean across proteins of medians of paired solve-time ratios; only fully observed three-arm cases",
  "resampling": "joint repetition indices across all three ratios, independently within each protein/lambda/thread case",
  "interval": "linear-interpolated 2.5/97.5 percentile interval; absent when draws=0 or no complete cases",
  "limitation": "repeat variability conditional on explicitly listed cases; incomplete subsets are not full-panel estimates; no population generalization or correction for systematic drift"
}
```
