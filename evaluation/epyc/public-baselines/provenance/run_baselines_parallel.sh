#!/usr/bin/env bash
# Native comparison as 8 parallel chains (user, 2026-09-25): 4 tools x 2 lambdas, each chain one
# complete runner invocation (feasibility, then 5 repetitions), pinned to its own CPU with memory
# bound to that CPU's node. At most 8 cores: 2 on node 0, 6 on node 1, one per L3 group.
# Waits for the j16 rerun launcher to finish first, then merges the chains and runs the analysis.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="$REPO/research/results/epyc-node1-2026-09"
TIMING="$OUT/timings"
TOOLS="$OUT/tools"
ROOT="$OUT/public-baselines"
BASELINE_TIMEOUT=240   # i9 cap 120 s x 2, the same factor as the timing-study caps

while pgrep -f 'rerun_j16_after_[p]arallel' > /dev/null; do sleep 30; done
echo "[$(date -u +%FT%TZ)] start native comparison, 8 parallel chains" | tee -a "$OUT/run_all.log"

chain() {  # chain <tool> <lambda> <cpu>
  local tool=$1 lam=$2 cpu=$3 node=$(( $3 / 16 ))
  local dir="$ROOT/chains/$tool-l$lam"
  mkdir -p "$dir"
  local common=(--panel "$TIMING/inputs/panel.json" --table "$TIMING/source/data/codon_usage_freq_table_human.csv"
    --sources "$OUT/baseline-sources" --lock "$dir/compute.lock" --tools "$tool" --lambdas "$lam"
    --sparsedesign-binary "$TIMING/bin/sparsedesign-square" --sparsedesign-source "$TIMING/source"
    --timeout "$BASELINE_TIMEOUT" --memory-gib 32 --cpu "$cpu")
  "$TOOLS/numabind" "$node" python3 "$TOOLS/publication_baselines.py" "${common[@]}" \
    --output "$dir/feasibility" --repetitions 1 --stop-after-cap >> "$dir/chain.log" 2>&1
  "$TOOLS/numabind" "$node" python3 "$TOOLS/publication_baselines.py" "${common[@]}" \
    --output "$dir/repeated" --repetitions 5 --feasibility "$dir/feasibility" >> "$dir/chain.log" 2>&1
  echo "[$(date -u +%FT%TZ)] chain $tool-l$lam done" | tee -a "$OUT/run_all.log"
}

# Node 1: one core per L3 group, then two second cores; node 0: two cores on separate L3 groups.
chain lineardesign  0 16 & chain lineardesign  4 20 &
chain linearcdsfold 0 24 & chain linearcdsfold 4 28 &
chain derna         0 17 & chain derna         4 21 &
chain sparsedesign  0 8  & chain sparsedesign  4 12 &
wait

python3 "$REPO/research/epyc-node1/merge_baselines.py" --root "$ROOT" | tee -a "$OUT/run_all.log"
exec "$REPO/research/epyc-node1/run_all.sh" analysis
