#!/usr/bin/env bash
# Re-run the i9-14900KF benchmarks on EPYC 7313 NUMA node 1 (CPUs 16-31, memory node 1).
# Prepared 2026-09-23; see AGENTS.md "EPYC node-1 replication". Resumable: rerun the same
# command after an interruption and completed records are verified and skipped.
#
# Usage (from the repo root, inside tmux/screen so a dropped SSH session does not kill it):
#   research/epyc-node1/run_all.sh [phase ...]
# Phases, in default order:
#   parallel baselines-feasibility baselines-repeated analysis
# 'parallel' (default since 2026-09-25) runs every remaining single-thread timing and follow-up task
# in parallel on up to 10 node-1 CPUs, scaled to host load (run_parallel.py). The multi-thread
# timing tasks were already complete; 'timings' and 'followup' still run them serially if named.
# 'timings' runs the four timing-study phases the paper reports (ablation profile scaling workstation),
# merged and
# ordered by descending thread count (run_by_threads.py); the single phases can still be named.
# The optional 'july' phase (July i9 lineardesign-clean runs) is not in the default list: those runs
# are not used by the paper in publication/. Name it explicitly to run it.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="$REPO/research/results/epyc-node1-2026-09"
TIMING="$OUT/timings"
FOLLOWUP="$OUT/followup-Q8VIM6"
TOOLS="$OUT/tools"
NUMABIND="$TOOLS/numabind"
LOG="$OUT/run_all.log"
BASELINE_TIMEOUT=240   # i9 cap 120 s x 2, the same factor as the timing-study caps

PHASES=("$@")
if [ ${#PHASES[@]} -eq 0 ]; then
  PHASES=(parallel baselines-feasibility baselines-repeated analysis)
fi

say() { echo "[$(date -u +%FT%TZ)] $*" | tee -a "$LOG"; }

# One launcher at a time.
exec 9>"$OUT/run_all.lock"
if ! flock -n 9; then echo "Another run_all.sh owns $OUT" >&2; exit 1; fi

[ "$(hostname -s)" = arrakis ] || { echo "Expected host arrakis, got $(hostname -s)" >&2; exit 1; }
for f in "$TIMING/epyc-node1-adaptation.json" "$FOLLOWUP/protocol.json" "$NUMABIND"; do
  [ -e "$f" ] || { echo "Missing $f; run the preparation steps in AGENTS.md first" >&2; exit 1; }
done

collector() {  # collector <study-dir> <phase>
  python3 "$1/source/publication_campaign_final.py" run --output "$1" --phase "$2" 2>&1 | tee -a "$1/$2-progress.log"
}

baselines() {  # baselines <feasibility|repeated>
  local extra=(--repetitions 1 --stop-after-cap)
  if [ "$1" = repeated ]; then extra=(--repetitions 5 --feasibility "$OUT/public-baselines/feasibility"); fi
  # numabind binds the runner (and so every native child) to node-1 memory and CPUs 16-31;
  # the runner then pins each native process to CPU 16.
  "$NUMABIND" 1 python3 "$TOOLS/publication_baselines.py" \
    --panel "$TIMING/inputs/panel.json" \
    --table "$TIMING/source/data/codon_usage_freq_table_human.csv" \
    --sources "$OUT/baseline-sources" \
    --output "$OUT/public-baselines/$1" \
    --lock "$OUT/compute.lock" \
    --tools lineardesign linearcdsfold derna sparsedesign \
    --sparsedesign-binary "$TIMING/bin/sparsedesign-square" \
    --sparsedesign-source "$TIMING/source" \
    --timeout "$BASELINE_TIMEOUT" --memory-gib 32 --cpu 16 "${extra[@]}" 2>&1 | tee -a "$OUT/baselines-$1.log"
}

for phase in "${PHASES[@]}"; do
  say "start $phase (load: $(cut -d' ' -f1-3 /proc/loadavg))"
  case "$phase" in
    timings) python3 "$REPO/research/epyc-node1/run_by_threads.py" --output "$TIMING" 2>&1 | tee -a "$TIMING/by-threads-progress.log" ;;
    pilot|ablation|profile|scaling|workstation) collector "$TIMING" "$phase" ;;
    parallel) python3 "$REPO/research/epyc-node1/run_parallel.py" --study "$TIMING" --study "$FOLLOWUP" 2>&1 | tee -a "$OUT/parallel-progress.log" ;;
    followup) python3 "$REPO/research/epyc-node1/run_by_threads.py" --output "$FOLLOWUP" --phases ablation 2>&1 | tee -a "$FOLLOWUP/by-threads-progress.log" ;;
    baselines-feasibility) baselines feasibility ;;
    baselines-repeated) baselines repeated ;;
    july) python3 "$REPO/research/epyc-node1/july.py" run 2>&1 | tee -a "$OUT/july.log" ;;
    analysis)
      # No --require-complete: the 24 pilots are deliberately not all run (not in the paper). Check the
      # analyzer's missing list: it should contain pilot tasks only.
      python3 "$TIMING/analyze_node1.py" --input "$TIMING" --output "$OUT/analysis/timings" \
        2>&1 | tee "$OUT/analysis-timings.log"
      python3 "$TOOLS/analyze_publication_baselines.py" --runs "$OUT/public-baselines/feasibility" "$OUT/public-baselines/repeated" \
        --panel "$TIMING/inputs/panel.json" --table "$TIMING/source/data/codon_usage_freq_table_human.csv" \
        --lock "$OUT/compute.lock" --output "$OUT/analysis/public-baselines" 2>&1 | tee "$OUT/analysis-baselines.log"
      ;;
    *) echo "Unknown phase $phase" >&2; exit 2 ;;
  esac
  say "done $phase"
done
say "all requested phases complete"
