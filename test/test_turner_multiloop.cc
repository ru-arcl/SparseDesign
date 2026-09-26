// Production-level regression whose returned optimum actually contains a Turner multiloop.
// The smaller exhaustive tests validate global objectives but happen to select only unpaired or
// single-stem structures; this case exercises sparse M2 through an enclosing closed pair.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "../src/codon_table.h"
#include "../src/dfa.h"
#include "../src/fold_turner.h"

// The vendored Vienna parameter table defines this mutable symbol. Changing it only inside this
// regression lets us exercise parameter sets whose per-unpaired multiloop term is nonzero without
// adding a test-only setter to the production energy API.
extern int ML_BASE37;

using namespace ldclean;

static double codon_penalty(const CodonTable& table, const std::string& rna) {
    double penalty = 0.0;
    for (size_t i = 0; i < rna.size(); i += 3)
        penalty -= std::log(table.w_of(rna.substr(i, 3)));
    return penalty;
}

// Return the largest number of unpaired residues after the last branch of any multiloop, or -1 if
// there is no multiloop. A positive result specifically exercises m1's direct-branch-plus-suffix
// shortcut, the case in which the dense reference must add the per-unpaired Turner term itself.
static int max_multiloop_trailing(const std::string& structure) {
    std::vector<int> stack;
    std::vector<int> close(structure.size(), -1);
    for (size_t i = 0; i < structure.size(); ++i) {
        if (structure[i] == '(') {
            stack.push_back(static_cast<int>(i));
        } else if (structure[i] == ')') {
            if (stack.empty()) return -1;
            close[stack.back()] = static_cast<int>(i);
            stack.pop_back();
        } else if (structure[i] != '.') {
            return -1;
        }
    }
    if (!stack.empty()) return -1;

    int max_trailing = -1;
    for (size_t open = 0; open < close.size(); ++open) {
        if (close[open] < 0) continue;
        int children = 0;
        int last_child_close = static_cast<int>(open);
        for (int cursor = static_cast<int>(open) + 1; cursor < close[open];) {
            if (structure[cursor] == '(') {
                ++children;
                last_child_close = close[cursor];
                cursor = close[cursor] + 1;
            } else {
                ++cursor;
            }
        }
        if (children >= 2)
            max_trailing = std::max(max_trailing, close[open] - last_child_close - 1);
    }
    return max_trailing;
}

static bool verify_case(const CodonTable& table, const std::string& protein, int ml_base,
                        bool require_trailing) {
    const double lambda = 400.0;  // CLI/paper lambda 4.
    ML_BASE37 = ml_base;
    DFA lattice(protein, table);

    TurnerSparseStats stats;
    TurnerDesign sparse = design_turner(lattice, table, lambda, 3, 1, &stats);
    double dense = turner_min_cost_dense_reference(lattice, lambda, 3);
    bool same_cost = std::fabs(sparse.cost - dense) < 1e-6;
    int trailing = max_multiloop_trailing(sparse.structure);
    bool has_multiloop = trailing >= 0;
    bool trailing_ok = !require_trailing || trailing > 0;
    bool exercised_candidates = stats.candidates > 0 && stats.direct_intervals > stats.candidates;
    bool shape_ok = sparse.mrna.size() == 3 * protein.size() &&
                    sparse.structure.size() == sparse.mrna.size();
    bool translates = shape_ok;
    if (translates)
        for (size_t i = 0; i < sparse.mrna.size(); i += 3)
            if (table.aa_of(sparse.mrna.substr(i, 3)) != protein[i / 3]) translates = false;
    double realized = shape_ok
        ? turner_eval(sparse.mrna, sparse.structure) + lambda * codon_penalty(table, sparse.mrna)
        : 1e100;
    bool backtrace_ok = std::fabs(sparse.cost - realized) < 1e-6;
    bool ok = same_cost && has_multiloop && trailing_ok && exercised_candidates && translates &&
              backtrace_ok;

    std::printf("Turner multiloop optimum (ML_BASE37=%d): dense=%.9f sparse=%.9f "
                "realized=%.9f Z=%zu direct=%zu trailing=%d translates=%s %s\n",
                ml_base, dense, sparse.cost, realized, stats.candidates, stats.direct_intervals,
                trailing, translates ? "yes" : "no", ok ? "OK" : "FAIL");
    return ok;
}

int main() {
    CodonTable table("data/codon_usage_freq_table_human.csv");
    int saved_ml_base = ML_BASE37;
    bool zero = verify_case(table, "WIPKKSAWCFEQQLCNLYSQLFMGV", 0, false);
    // This nearby sequence has an optimal multiloop with a two-base trailing unpaired suffix.
    // A nonzero value makes omission of that suffix penalty observable in the dense oracle.
    bool nonzero = verify_case(table, "WIPKKSAWCFEQQLCNLYSQGFMGV", 37, true);
    ML_BASE37 = saved_ml_base;
    return zero && nonzero ? 0 : 1;
}
