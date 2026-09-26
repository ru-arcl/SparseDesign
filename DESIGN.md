# SparseDesign algorithm and implementation

SparseDesign implements exact weighted-codon-automaton design under the Turner-d0 model.
The implementation lineage is recorded in [PROVENANCE.md](PROVENANCE.md); original code uses
[Rutgers Non-commercial Research License (RU-NCRL)](LICENSE.md), with separate third-party notices.

The published LinearDesign paper and supplement supply the weighted-DFA formulation
(Zhang et al., Nature 621:396–403, 2023; DOI 10.1038/s41586-023-06127-z). Established sparse
RNA-folding recurrences motivate candidate normalization. The accompanying SparseDesign
manuscript gives the complete multiloop recurrence and relative-exactness proof, including
concrete endpoint-state joins and the scalar parent interface.

Here, exactness means preservation of the specified dense recurrence in real arithmetic.
The implementation uses double precision; its regression and reconstruction checks assess
the numerical implementation, rather than proving bit-identical evaluation of every recurrence.
This document describes the released implementation. References to the LinearDesign supplement
below provide background; the SparseDesign manuscript states the sparsification theorem.

## Problem and objective

Given a protein `p = p₀ … p_{|p|−1}`, design an mRNA CDS `r` (a synonymous coding sequence)
that minimizes the joint objective:

```
MFECAI_λ(r) = MFE(r) − (|r|/3)·λ·log CAI(r)
            = MFE(r) − λ · Σᵢ log w(codon(r, i))          (expanded; |r|/3 = #codons)
r*_λ(p) = argmin over r encoding p of MFECAI_λ(r)
```

- `MFE(r) = min over allowed secondary structures s of ΔG°(r, s)` — RNA folding free energy;
  structures are pseudoknot-free, have hairpin loops of at least three nucleotides, and have
  at most 30 unpaired nucleotides in each bulge/internal loop.
- `w(c)` = **relative adaptiveness** of codon `c` = freq(c) / max freq over synonymous codons
  of the same amino acid (0 < w ≤ 1). A zero frequency, including an all-zero synonymous
  family, is assigned `w = 1e-9`; zero does not exclude a codon.
- `CAI(r) = (∏ᵢ w(codon_i))^(1/#codons)` — geometric mean (Sharp & Li 1987).
- `λ ≥ 0` trades MFE against codon optimality. `λ = 0` gives pure MFE; larger λ gives
  greater weight to CAI, and can select the same optimum over an interval of λ values.
- The implementation designs exactly the supplied residues. It does not append a stop codon;
  include `*` explicitly in the protein when a stop is part of the desired CDS.

Because the objective adds **−λ·log w(c)** per codon, codon optimality is *decomposable per
codon* and folds directly into the folding DP as additive weights.

## Energy model

- **Production model:** Turner 2004 nearest-neighbor energies with dangling ends disabled
  (**dangle model 0**, Vienna `-d0`) and special tri-, tetra- and hexaloops enabled.
  The energy functions follow published loop-energy formulas; the structure evaluator is
  a separate traversal sharing those functions with the optimizer.
- Energies are integers in units of 0.01 kcal/mol (Vienna convention), at 37°C. The CAI
  term is scaled to the same units but remains double precision; CLI `λ` is multiplied by 100.
- **Simplified model (validation only):** the LinearDesign supplement's pseudocode (Figs. 2–4) uses a
  Nussinov–Jacobson model (per-base-pair: ΔG(C,G)=−3, ΔG(A,U)=−2, ΔG(G,U)=−1). We implement
  this model in `src/fold_simple.*` to check the lattice, DFA, DP and backtrace machinery
  against exhaustive enumeration on short proteins. It is not a CLI design mode.

## Weighted codon DFA

- DFA `D = (Q, Σ, δ, q₀, F)`, `Σ = {A,C,G,U}`, edges labeled by a nucleotide, start `q₀=(0,0)`.
- Per amino acid, a small prefix-trie DFA whose start→end paths spell exactly its synonymous
  codons. Shared prefixes naturally accommodate the 6-codon Leu/Arg/Ser families.
- Concatenate per-protein: `D(p) = D(p₀) ∘ D(p₁) ∘ … ∘ D(p_{|p|−1})`, stitching each DFA's end
  state to the next's start. RNA length is `n = 3|p|`; lattice width is bounded by the genetic
  code for the unmodified codon DFA.
- **CAI weights on edges:** distribute `−log w(c)` across each codon's edges so every codon
  path sums to `−log w(c)`. We place the full `−log w(c)` on the codon's **third-nucleotide
  edge** (first/second edges weight 0); since synonymous codons within a prefix group differ
  at the third position, each third-edge maps to a unique codon.
- `out_edges(q)`, `in_edges(q)` return `(neighbor, nucleotide, weight)` triples.

Note: DFA *minimization* affects resource use, not the optimal value — the optimum is invariant
to how the synonymous set is encoded, as long as it is encoded exactly with correct weights.

Repeated `--forbid-motif` options compose the codon DFA with a deterministic motif matcher.
Its accepted paths exclude the motifs throughout the CDS, including across codon boundaries.
The product preserves codon weights but can increase lattice width; bounds stated in RNA
length assume bounded width, including a fixed constraint automaton.

## Lattice parsing and folding DP

RNA folding = CKY parsing of a structure grammar; mRNA design = **lattice parsing** = the same
DP generalized from a single sequence to the DFA (Bar-Hillel intersection of the folding SCFG
with the weighted DFA). DP items are pairs of DFA nodes `(q_i, q_j)` instead of positions
`(i, j)`. The exact solver remains worst-case cubic in sequence length and quadratic in memory.

The LinearDesign supplement describes two formulations:

- **Bottom-up** (Fig. 2): CKY over spans, nonterminals S (any structure) and P (a paired
  region), rules `S→S N` (append unpaired), `P→a S b` (pair), `S→S P` (bifurcation).
- **Left-to-right** (Fig. 4): recurrences explored incrementally, with **beam search**
  (BeamPrune, Fig. 3) available for approximate search. SparseDesign uses a gap-ordered
  wavefront and retains the exact state space.

For the **Turner** model the single `ΔG(a,b)` pairing score is replaced by the loop-type
decomposition (hairpin / stack / bulge / internal / multiloop / external) using the Turner
energy functions, with dangle model 0 — i.e. the standard ViennaRNA/LinearFold V/M/M2/External
state decomposition, generalized to lattice nodes. CAI edge weights `λ·w_edge` are added as
each edge is traversed.

Backtrace recovers an optimal `(sequence, dot-bracket structure)`.

## Exact sparse multiloop recurrence

The production solver's cubic hot spot is the Turner multiloop decomposition. For a lattice-state
interval `(a,b)`, let `B(a,b)` be the cheapest direct closed branch including its dangle-0 branch
penalty, and let `P(a,b)` be the best partitionable or endpoint-unpaired fragment. A direct branch
is retained as a candidate exactly when `B(a,b) < P(a,b)`. The `M2` recurrence scans only retained
rightmost candidates ending at `b`, joining at the same concrete DFA state. A non-candidate
direct branch remains available to `M1`; only its use as a later split candidate is pruned.

Reassociating a non-candidate through its no-more-expensive `P` realization proves value exactness
under the dangle-0 interface-separability condition. With `N = |Q|` lattice states, edge set `E`
and `Z` retained candidates, multiloop work is `O(N² + N|E| + NZ)`, or `O(n²+nZ)` under
bounded lattice width. `Z` can be quadratic, so the complete solver's worst case
remains `O(n³)`. Candidate storage is `O(|Q|+Z)`, while the complete solver retains its `O(|Q|²)`
tables. On the measured 64-bit build, the candidate index uses 16 bytes per allocated candidate plus
a 24-byte vector header per right endpoint; the unchanged scalar/vector/flag cell objects alone use
about 74 bytes per ordered lattice-state pair before their heap payload. See the
accompanying manuscript for the theorem and full model/interface assumptions.

The scalar branch collapse is specific to the implemented Turner dangle model 0: a parent
multiloop observes a branch's pair type only through the already-included branch penalty. Dangles
or coaxial stacking would expose additional interface context and require richer candidates.

## Controlled performance attribution

`turner_min_cost_dense_right_normal` shares the production scalar right-normal recurrence,
collapsed branch cost, endpoint lists and iteration order. A compile-time retention switch
keeps every finite direct interval instead of applying `B < P`. This isolates the effect
of pruning. The older `turner_min_cost_dense_wavefront` is a left-normal pair-indexed
recurrence, so comparison with it includes changes to orientation, scalar collapse and
branch representation. It must not be described as pruning alone.

`turner_min_cost_profiled` uses separately instantiated instrumented kernels. It records
memo-construction time, three phase wall times (including their OpenMP barriers), and
split-list visits. `split_eligible` counts entries surviving endpoint-position guards,
before testing the left fragment's feasibility. Per-row counters have a unique writer
within a wavefront, avoiding atomics. Performance timing should use ordinary kernels.
Inner-loop counters can substantially perturb the profiled kernels, including their compiler
and memory behavior. Counts measure scanned work; phase times describe the instrumented
executions and are not a quantitative phase breakdown of ordinary-kernel elapsed time.

## Output

The `sparsedesign` executable reports the designed mRNA CDS, its dot-bracket structure,
`MFE` in kcal/mol (= integer energy / 100), and `CAI = exp((1/#codons)·Σ log w)`.
The legacy executable name `lineardesign-clean` is retained for script compatibility.
Main CLI options are `-l/--lambda`, `-c/--codonusage`, and
`-j/--threads`. Repeated `--forbid-motif` options compose the codon lattice with a finite-state
motif filter, and `--sparse-stats` emits candidate counts after an optional read-only post-pass.
The standard summary rounds MFE to two decimals and CAI to three decimals.
`--diagnostics` prints the full-precision DP objective and reconstructed objective to stderr.
Both objectives use 0.01 kcal/mol units; the separately reported codon penalty is dimensionless.
Reconstruction uses a separate structure traversal that shares the optimizer's energy functions.
Every design checks that the realized sequence follows the DFA and that its structure energy
plus DFA path penalty agrees within `1e-4 + 1e-10*max(abs(cost),abs(realized_cost))` internal units.
The CLI rejects a second FASTA record instead of concatenating unrelated inputs.
`--evaluate-mrna` validates and folds one supplied synonymous sequence using a fixed-path DFA;
it does not optimize codons and cannot be combined with nonzero lambda or motif constraints.

## Memo-layout builds

The production Makefile defines `LDCLEAN_SQUARE_MEMOS` for the square memo arena.
This was the measured winner at both one and sixteen threads on the recorded 1,006-aa Q6UXY8 screen. The
explicit `make lowmem-packed` target omits that definition and builds the forward-only,
position-major arena; it preserves the exact recurrence and output contract while reducing peak
RSS. Its runtime effect depends on input and hardware; the current Dp427c study records
a different layout tradeoff from the older screen. `make test-layouts` runs the exact sparse, motif, multiloop, wavefront, and
CLI checks under both layouts, and the sanitizer targets cover each build separately. Neither
layout changes the asymptotic quadratic state bound.

## Implemented components

1. Genetic code, codon-usage parsing, `w(c)`, and CAI.
2. Weighted codon DFA construction plus exact finite-state motif composition.
3. Simplified Nussinov lattice DP and brute-force validation.
4. Turner dangle-0 lattice DP, exact backtrace, compact forward state, and separate structure evaluation.
5. Exact candidate-sparse multiloops and a deterministic gap-ordered OpenMP wavefront.
6. Dense, exhaustive-CDS, exact-integer kernel, constrained-language, sanitizer, and long-sequence
   regression tests.

## Verification

Historical external validation used recorded outputs and black-box runs of the original
LinearDesign implementation. Its restricted binaries are not a build or test dependency of
this release. The comparison target is the optimal joint objective under aligned settings.
An equal weighted objective need not imply an identical sequence, structure, or separate
MFE and CAI values when optima tie.

The current suite additionally enumerates 1,552 synonymous CDSs across five proteins and compares
fixed-sequence dense folding, a dense weighted-lattice oracle, and the production sparse solver at
four objective weights. An exact-integer abstract kernel checks every `M1/M2` cell on 500 generated
branching/merging DAGs, motif-product tests verify the accepted language explicitly, and a selected
production optimum exercises a true multiloop under both zero and nonzero per-unpaired multiloop
penalties.

The release also compares every M1/M2 and closed-pair cell of the three actual Turner
recurrences on 72 configurations (181,344 intervals), using motif products, three
lambda values, nonzero/zero multiloop-unpaired penalties, and seeded random proteins.
All-cell accessors are compiled only into this test; no debug API is exported.
