# Commodity-PC Dp427c demonstration

This compact dataset retains the ten i9-14900KF Dp427c measurements cited as a
separate commodity-PC example in the SparseDesign paper. The primary benchmark
platform is [AMD EPYC](../epyc/README.md); results from the two platforms are kept
separate.

Dp427c (NP_000100.3) has 3,677 amino acids, giving an 11,031-nucleotide coding
sequence. All runs use lambda 0, the human codon table, and 16 threads. The host
has approximately 64 GiB physical RAM and a hybrid P/E-core topology; this
16-thread configuration includes both types. GCC 14.3.0 compiled the measured
executables with `-O3 -flto` and OpenMP. Each process had a 32-GiB virtual-address
limit and an 1,800-second wall-time cap. The address-space limit does not describe
the host's physical memory.

| Memo layout | Successful runs | Median wall time | Median peak RSS |
|---|---:|---:|---:|
| Square | 5 | 128.26 s | 20.10 GiB |
| Packed | 5 | 126.42 s | 14.43 GiB |

Wall time and peak RSS cover the complete command and come from GNU time.
They are not kernel-only measurements. The [CSV](analysis/workstation.csv) gives
the full per-layout summaries, and [JSON](analysis/summary.json) also preserves
the five paired-layout comparisons.

## Contents and provenance

The raw `record.json`, stdout, stderr, and GNU-time files for all ten runs are
under `timings/runs/`. These files, the Dp427c FASTA input, frozen source,
build record, CPU metadata, and original protocols were copied byte-for-byte
from the September 2026 i9 timing dataset. `SUBSET_ORIGIN.json` maps each copied
file to its original path and SHA-256. The JSON summary contains only the
`workstation` and `workstation_paired` sections extracted from the original
analysis; the CSV is an unchanged copy.

The original protocol and panel metadata include registrations for other phases
so their recorded hashes remain independently checkable. Only the Dp427c input
and run outputs are retained here. The source snapshots and absolute command
paths preserve collection provenance; they are not a launcher for the omitted
studies. Optimized executables are not distributed.

## Verify and regenerate the summary

From the standalone repository root, with Python 3.10 or later:

```bash
python3 -B evaluation/commodity-pc/verify.py
python3 -B evaluation/commodity-pc/verify.py --output /tmp/commodity-pc-summary.json
```

The verifier checks the complete dataset manifest, original copied-file hashes,
source and input identities, task registrations, successful outcomes, output
hashes, and agreement between raw GNU-time fields and the records. It then
recomputes both summary sections and checks the recorded JSON and CSV. It runs
no optimizer and requires no third-party Python packages. The optional output
must be outside this dataset to preserve its checksums.
