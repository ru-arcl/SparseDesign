// fold_simple.h — simplified Nussinov-Jacobson lattice DP (validation milestone).
// Clean-room: the simplified energy model and lattice-parsing structure are from the
// LinearDesign paper/supplement (Sup Figs 2-4): per-base-pair scores ΔG(C,G)=-3,
// ΔG(A,U)=-2, ΔG(G,U)=-1, joint with the CAI cost λ·Σ(-log w). This validates the
// lattice DP + backtrace machinery against brute force before the full Turner model.
#ifndef LDCLEAN_FOLD_SIMPLE_H
#define LDCLEAN_FOLD_SIMPLE_H

#include <string>

#include "dfa.h"

namespace ldclean {

struct DesignResult {
    std::string mrna;       // designed coding sequence (with stop)
    std::string structure;  // dot-bracket
    double cost;            // total objective: folding energy + λ·Σ(-log w)
    double fold_energy;     // simplified Nussinov folding energy component
    double cai_cost;        // λ·Σ(-log w) component
};

// Exact (no beam) simplified-model design over the DFA. min_loop = min unpaired in a hairpin.
DesignResult design_simple(const DFA& dfa, double lambda, int min_loop = 3);

// Simplified Nussinov base-pair score for two nucleotides (1..4), or 0 if disallowed.
int simple_pair_score(int nuc_a, int nuc_b);
bool simple_pairable(int nuc_a, int nuc_b);

}  // namespace ldclean

#endif  // LDCLEAN_FOLD_SIMPLE_H
