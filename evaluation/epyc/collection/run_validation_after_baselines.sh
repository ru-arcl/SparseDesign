#!/usr/bin/env bash
# 2026-09-25 (user): rerun the i9 validation study on EPYC node 1, one core (CPU 16), node-1 memory.
# It is the source of the paper's "about 15.6 s / 109.5 MiB" finite-state-constraint figure. Starts once
# the parallel native comparison has been merged, so it runs alongside the analysis (2 cores in use).
# Frozen script/driver/solver sources match the i9 build hashes; the compiler here is GCC 11.4 (i9: 14.3).
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="$REPO/research/results/epyc-node1-2026-09"
V="$OUT/validation"
until [ -f "$OUT/public-baselines/repeated/metadata.json" ]; do sleep 60; done
echo "[$(date -u +%FT%TZ)] start validation (CPU 16)" | tee -a "$OUT/run_all.log"
"$OUT/tools/numabind" 1 taskset -c 16 python3 "$REPO/publication/code/evaluation/validation/source/publication_validation.py" \
  --root "$REPO/publication/code" --source "$OUT/timings/source" --out "$V" --build --lock "$V/compute.lock" \
  >> "$V/validation.log" 2>&1 && status=ok || status="failed ($?)"
echo "[$(date -u +%FT%TZ)] validation $status" | tee -a "$OUT/run_all.log"
