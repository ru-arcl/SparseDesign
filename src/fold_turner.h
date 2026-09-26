// fold_turner.h — SparseDesign's Turner-d0 lattice dynamic program.
// Clean-room: lattice-parsing structure from the LinearDesign paper (Sup Figs 2-4) extended
// from the simplified model to the Turner nearest-neighbor loop decomposition (hairpin /
// stack / bulge / internal / multiloop / external), dangle model 0, using the open energy
// module (energy.h). Joint objective: folding energy + λ·Σ(-log w).
#ifndef LDCLEAN_FOLD_TURNER_H
#define LDCLEAN_FOLD_TURNER_H

#include <cstddef>
#include <string>

#include "codon_table.h"
#include "dfa.h"

namespace ldclean {

struct TurnerDesign {
    std::string mrna;       // designed coding sequence
    std::string structure;  // dot-bracket MFE structure
    double cost;            // total objective (folding energy in 0.01 kcal/mol + λ·CAI cost)
    double mfe_kcal;        // folding free energy in kcal/mol
    double cai;             // CAI of the designed sequence
    double realized_cost;   // independently rescored structure + lambda * realized DFA path cost
    double codon_penalty;   // sum of -log w(codon), using the supplied table (also for fixed RNA)
};

// Optional publication/diagnostic counters for the exact sparse multiloop recurrence. Collecting
// these counters performs an O(num_nodes^2) read-only post-pass; normal solves pay no scan cost.
struct TurnerSparseStats {
    int rna_length = 0;
    std::size_t lattice_nodes = 0;
    std::size_t interval_cells = 0;       // layered state pairs filled by the wavefront
    std::size_t direct_intervals = 0;     // cells with at least one feasible direct closed branch
    std::size_t candidates = 0;           // retained exact multiloop candidates Z
    std::size_t candidate_capacity = 0;   // allocated MultiCandidate slots
    std::size_t candidate_storage_bytes = 0;  // vector headers + allocated candidate slots
    std::size_t max_candidates_per_end = 0;
};

enum class TurnerKernel { Sparse, DenseRightNormal, DenseWavefront };

// Separate instrumented instantiations count split-list visits and wall time at phase barriers.
// Production solves do not execute counters or clock reads. Timings include OpenMP overhead;
// construction_seconds measures memo allocation/initialization but excludes DFA construction.
struct TurnerProfile {
    double construction_seconds = 0.0;
    double closed_seconds = 0.0;
    double multiloop_seconds = 0.0;
    double external_seconds = 0.0;
    std::size_t split_visits = 0;
    std::size_t split_eligible = 0;
};

// Minimum joint cost over the design space (0.01 kcal/mol units + λ·CAI). At λ=0 this is the
// pure minimum folding energy (MFE*100). min_loop = minimum hairpin loop size (Turner = 3).
// Throws std::invalid_argument when lambda is non-finite/negative or min_loop is negative.
double turner_min_cost(const DFA& dfa, double lambda, int min_loop = 3, int threads = 1,
                       TurnerSparseStats* sparse_stats = nullptr);

// Cubic dense multiloop recurrence retained as a validation oracle. This is intentionally
// single-threaded and is not used by the production solver.
double turner_min_cost_dense_reference(const DFA& dfa, double lambda, int min_loop = 3);

// Validation/ablation path: evaluate the same cubic dense multiloop recurrence with the current
// compact bottom-up wavefront, eager contiguous branch rows, and fused m1/m2 kernel. This measures
// the whole multiloop replacement, including orientation and pair-type collapse. Production calls do
// not use this path; its eager branch rows require O(B) extra memory, where B is the number of
// feasible direct branch/pair-type entries (O(num_nodes^2) in the worst case).
double turner_min_cost_dense_wavefront(const DFA& dfa, double lambda, int min_loop = 3,
                                       int threads = 1);

// Controlled ablation: the scalar right-normal recurrence, endpoint lists and pair-type collapse
// are identical to production, but EVERY finite direct interval is retained. Only the retention
// predicate differs. A sparse_stats argument reports the retained dense list, not a pruned Z.
double turner_min_cost_dense_right_normal(const DFA& dfa, double lambda, int min_loop = 3,
                                          int threads = 1, TurnerSparseStats* sparse_stats = nullptr);

double turner_min_cost_profiled(const DFA& dfa, double lambda, TurnerKernel kernel,
                                TurnerProfile& profile, int min_loop = 3, int threads = 1,
                                TurnerSparseStats* sparse_stats = nullptr);

// Full design: optimal CDS + dot-bracket structure + MFE + CAI, via backtrace. The exact sparse
// bottom-up fill is serial at threads=1 and parallel across each span diagonal at threads>1.
// Throws std::overflow_error instead of attempting a backtrace when the objective overflows.
TurnerDesign design_turner(const DFA& dfa, const CodonTable& table, double lambda,
                           int min_loop = 3, int threads = 1,
                           TurnerSparseStats* sparse_stats = nullptr);

// Runtime team size observed for the solver's capped `requested` value. The policy ceiling is all
// cores on a single-socket host and one node's cores on a multi-socket host (spanning sockets
// regresses this bandwidth-sensitive DP); the observed team can be smaller in a non-OpenMP build
// or under runtime limits/dynamic adjustment. Lets callers report the effective team honestly.
int turner_effective_threads(int requested);

// Independent evaluator: Turner-d0 free energy (0.01 kcal/mol) of a given sequence + structure.
// Used to verify that a designed (sequence, structure) actually scores the reported MFE.
int turner_eval(const std::string& seq, const std::string& structure);

}  // namespace ldclean

#endif  // LDCLEAN_FOLD_TURNER_H
