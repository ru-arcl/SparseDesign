# EPYC native public-software comparison

These are the primary native-software measurements for the manuscript: 96 feasibility
conditions and 210 fresh repetitions on the dual AMD EPYC 7313 server `arrakis` (1 TiB
RAM), built with GCC 11.4.0. Every native process used one thread, a 240-second wall cap
and a 32-GiB virtual-address-space cap. The superseded i9 native comparison is
omitted from this code release.

All 210 fresh repetitions succeeded. Feasibility contributed 42 successes, 8 timeouts
and 46 prescribed skips. The 252 successful executions independently validate and yield
36 distinct RNA/structure pairs (30 RNA sequences). Source-audited DERNA/SparseDesign
ratios at lambda 0 are 5.08, 4.71, 5.51 and 6.04 for 82, 95, 255 and 310 residues.
Only these four conditions meet all matched-comparison gates; native models otherwise
differ. See `protocol.json`, the recorded source identities and the analysis below.

`feasibility/` and `repeated/` preserve all retained raw records and native output bytes,
including task manifests, complete per-chain metadata, executable/source digests and
resource limits. `inputs/` contains the approved protein sequences and codon table.
`parameter-audit.json` binds the thermodynamic audit to the actual measured source
hashes. `provenance/` retains build evidence and the launcher/merge scripts.
Restricted upstream sources, executables and libraries are excluded.

Benchmarks ran on a shared server whose load varied during the measurements. Native
records retain initial load averages; they do not contain before/after CPU utilization
snapshots. The recorded metadata preserve the actual execution details.

## Verify and reproduce the analysis

From this directory, verify the complete file set and its recorded bytes:

```bash
python3 - <<'PY'
import hashlib, json
from pathlib import Path
root = Path('.')
expected = json.loads((root / 'sha256.json').read_text())
actual = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
          for p in root.rglob('*') if p.is_file()
          and p.name != 'sha256.json' and '__pycache__' not in p.parts}
assert actual == expected
print('Verified', len(actual), 'files')
PY
```

The analyzer requires ViennaRNA 2.7.2. It verifies input, source, executable, task and
raw-output digests, then independently translates and rescores retained outputs.
It does not run new design solves. Use an output directory outside this frozen export:

```bash
mkdir -p /tmp/sparsedesign-epyc-native-analysis
cp analysis/rescore-cache.json /tmp/sparsedesign-epyc-native-analysis/
python3 source/analyze_publication_baselines.py \
  --runs feasibility repeated --panel inputs/panel.json \
  --table inputs/human-codon-table.csv --model-audit parameter-audit.json \
  --lock /tmp/sparsedesign-epyc-native-analysis.lock \
  --output /tmp/sparsedesign-epyc-native-analysis
MPLCONFIGDIR=/tmp/sparsedesign-matplotlib python3 source/plot_publication_baselines.py \
  --analysis /tmp/sparsedesign-epyc-native-analysis \
  --output /tmp/sparsedesign-epyc-native-analysis --timeout 240
```

The cache contains independently computed scores under the recorded common model;
omit the cache copy to recompute these scores and fixed-sequence folds. File paths in
regenerated analysis may reflect the new location. Plot metadata may change with the
plotting environment; the measured values and 240-second timeout coordinates must agree.
The EPYC plotter changes the frozen i9 plotter only to derive caption text from the
`--timeout` value; that adaptation is recorded in `provenance/plot-adaptation.json`.

To repeat the source-parameter audit, independently obtain DERNA at the pinned revision
in `protocol.json`, then run `source/audit_publication_baseline_models.py --sources
<upstream-checkouts> --sparsedesign-source ../timings/source --output <new-audit.json>`.
The source hashes must match those in the measured metadata before using the audit.

From the publication `paper/` directory, import this evidence and regenerate numerical
claims with `python3 tools/import_epyc_native_evidence.py`. This verifies the complete
EPYC manifest before updating the local native data, figure, table and claim macros.
