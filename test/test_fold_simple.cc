// test_fold_simple.cc — brute-force validation of the simplified-Nussinov lattice DP.
// For short proteins we enumerate the ENTIRE synonymous design space, fold each sequence
// with a standard single-sequence Nussinov DP (same ΔG and min-loop), add the CAI cost,
// and confirm the lattice DP finds the same global-optimum objective value.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "../src/codon_table.h"
#include "../src/dfa.h"
#include "../src/fold_simple.h"

using namespace ldclean;

// Standard single-sequence Nussinov minimum folding energy (simplified ΔG, min-loop).
static double nussinov(const std::string& s, int min_loop) {
    int n = (int)s.size();
    std::vector<int> nc(n);
    for (int i = 0; i < n; ++i) nc[i] = nuc_code(s[i]);
    std::vector<std::vector<double>> dp(n + 1, std::vector<double>(n + 1, 0.0));
    for (int span = 1; span <= n; ++span) {
        for (int i = 0; i + span - 1 < n; ++i) {
            int j = i + span - 1;
            double best = dp[i + 1][j];  // i unpaired
            for (int k = i + min_loop + 1; k <= j; ++k) {
                if (!simple_pairable(nc[i], nc[k])) continue;
                double inner = (i + 1 <= k - 1) ? dp[i + 1][k - 1] : 0.0;
                double right = (k + 1 <= j) ? dp[k + 1][j] : 0.0;
                best = std::min(best, simple_pair_score(nc[i], nc[k]) + inner + right);
            }
            dp[i][j] = best;
        }
    }
    return dp[0][n - 1];
}

// Enumerate the synonymous design space; return min objective and the argmin sequence.
static void brute(const CodonTable& t, const std::string& residues, size_t k,
                  std::string& acc, double lambda, int min_loop,
                  double& best, std::string& best_seq) {
    if (k == residues.size()) {
        double cai_cost = 0.0;
        for (size_t i = 0; i < acc.size(); i += 3)
            cai_cost += -std::log(t.w_of(acc.substr(i, 3)));
        double obj = nussinov(acc, min_loop) + lambda * cai_cost;
        if (obj < best - 1e-12) { best = obj; best_seq = acc; }
        return;
    }
    for (const auto& opt : t.options_for(residues[k])) {
        acc += opt.codon;
        brute(t, residues, k + 1, acc, lambda, min_loop, best, best_seq);
        acc.erase(acc.size() - opt.codon.size());
    }
}

int main() {
    CodonTable t("data/codon_usage_freq_table_human.csv");
    const int min_loop = 3;
    int failures = 0;

    std::vector<std::string> proteins = {"MF", "MLF", "MWY", "MFLK", "MDEYW", "MKLVFR"};
    std::vector<double> lambdas = {0.0, 0.5, 1.0, 3.0};

    for (const auto& prot : proteins) {
        DFA dfa(prot, t);
        std::string residues = prot;  // no auto-stop
        double space = 1.0;
        for (char aa : residues) space *= t.options_for(aa).size();
        for (double lam : lambdas) {
            DesignResult r = design_simple(dfa, lam, min_loop);

            double bf_best = std::numeric_limits<double>::infinity();
            std::string bf_seq, acc;
            brute(t, residues, 0, acc, lam, min_loop, bf_best, bf_seq);

            // Self-consistency: the DP sequence's own Nussinov MFE + CAI must equal r.cost.
            double cai_cost = 0.0;
            for (size_t i = 0; i < r.mrna.size(); i += 3)
                cai_cost += -std::log(t.w_of(r.mrna.substr(i, 3)));
            double self = nussinov(r.mrna, min_loop) + lam * cai_cost;

            bool match = std::fabs(r.cost - bf_best) < 1e-6;
            bool selfok = std::fabs(r.cost - self) < 1e-6;
            bool transok = true;
            for (size_t i = 0; i < r.mrna.size(); i += 3)
                if (t.aa_of(r.mrna.substr(i, 3)) != residues[i / 3]) transok = false;

            printf("%-7s λ=%.1f space=%-6.0f  DP=%.4f  brute=%.4f  %s  self=%s  prot=%s\n",
                   prot.c_str(), lam, space, r.cost, bf_best,
                   match ? "MATCH" : "MISMATCH", selfok ? "ok" : "BAD",
                   transok ? "ok" : "BADTRANS");
            if (!match || !selfok || !transok) ++failures;
        }
    }

    // Show one example design.
    DFA dfa("MFLK", t);
    DesignResult r = design_simple(dfa, 1.0, min_loop);
    printf("\nexample MFLK λ=1.0:\n  %s\n  %s\n  cost=%.3f (fold=%.3f, cai_cost=%.3f) CAI=%.3f\n",
           r.mrna.c_str(), r.structure.c_str(), r.cost, r.fold_energy, r.cai_cost,
           t.cai(r.mrna));

    printf("\n%s\n", failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED");
    return failures == 0 ? 0 : 1;
}
