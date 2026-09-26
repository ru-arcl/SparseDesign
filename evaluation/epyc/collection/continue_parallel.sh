#!/usr/bin/env bash
# Replace the parallel scheduler without stopping runs in progress (2026-09-25, 8-core cap).
# The new scheduler adopts running workers, then the launcher continues with the remaining phases
# once every worker has exited (the workers hold the launcher lock inherited from the old launcher).
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="$REPO/research/results/epyc-node1-2026-09"
python3 "$REPO/research/epyc-node1/run_parallel.py" --study "$OUT/timings" --study "$OUT/followup-Q8VIM6" \
  >> "$OUT/parallel-progress.log" 2>&1
until flock -n "$OUT/run_all.lock" true; do sleep 10; done
exec "$REPO/research/epyc-node1/run_all.sh" baselines-feasibility baselines-repeated analysis
