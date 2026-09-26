#!/usr/bin/env bash
# 2026-09-25: wait for the last parallel worker to exit, then rerun the archived j16 tasks (thread order).
# Stops afterwards: the native comparison is started separately once the user decides how to run it.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
while pgrep -f 'run_parallel.py --[w]orker' > /dev/null; do sleep 10; done
exec "$REPO/research/epyc-node1/run_all.sh" timings
