// Compare every production M1/M2 and closed-pair cell against the original dense oracle.
// Include the implementation to expose test-only accessors without adding a public debug API.
#define LDCLEAN_TEST_CELL_PARITY
#include "../src/fold_turner.cc"
#include <cstdio>
#include <random>

extern int ML_BASE37;

using namespace ldclean;

static bool same(double a, double b) {
    return a == b || (std::isfinite(a) && std::isfinite(b) &&
                     std::fabs(a - b) <= 1e-7 + 1e-12 * std::max(std::fabs(a), std::fabs(b)));
}

int main() {
    CodonTable table("data/codon_usage_freq_table_human.csv");
    std::vector<std::string> proteins{"SSSS", "MSWFWFM", "WIPKKSAWCFEQQLCNLYSQLFMGV",
                                      "WIPKKSAWCFEQQLCNLYSQGFMGV"};
    std::mt19937 random(20260906);
    const std::string amino_acids = "ACDEFGHIKLMNPQRSTVWY";
    for (int k = 0; k < 8; ++k) {
        std::string protein;
        for (int j = 0; j < 6 + k; ++j) protein += amino_acids[random() % amino_acids.size()];
        proteins.push_back(protein);
    }
    std::size_t cells = 0;
    int cases = 0;
    int failures = 0;
    const int saved = ML_BASE37;
    for (const auto& protein : proteins) {
        for (double lambda : {0.0, 40.0, 400.0}) {
            for (int ml_base : {0, 37}) {
                ML_BASE37 = ml_base;
                // Every case exercises a product lattice; motifs may span codon boundaries.
                DFA dfa(protein, table, {"UCUUCU", "AGCAGC"});
                TurnerSolver sparse(dfa, lambda, 3, 1);
                TurnerSolver right(dfa, lambda, 3, 1);
                TurnerSolver dense(dfa, lambda, 3, 1);
                const double expected = dense.solve_dense_reference();
                if (!same(expected, sparse.solve()) || !same(expected, right.solve_dense_right_normal()))
                    ++failures;
                for (int p = 0; p < dfa.rna_length(); ++p)
                    for (int q = p + 1; q <= dfa.rna_length(); ++q)
                        for (int a : dfa.nodes_at(p))
                            for (int b : dfa.nodes_at(q)) {
                                bool ok = same(dense.test_m1(a, b), sparse.test_m1(a, b)) &&
                                          same(dense.test_m1(a, b), right.test_m1(a, b)) &&
                                          same(dense.test_m2(a, b), sparse.test_m2(a, b)) &&
                                          same(dense.test_m2(a, b), right.test_m2(a, b));
                                const auto& dv = dense.test_closed(a, b);
                                const auto& sv = sparse.test_closed(a, b);
                                const auto& rv = right.test_closed(a, b);
                                ok = ok && dv.size() == sv.size() && dv.size() == rv.size();
                                for (std::size_t k = 0; ok && k < dv.size(); ++k)
                                    ok = dv[k].nt5 == sv[k].nt5 && dv[k].nt3 == sv[k].nt3 &&
                                         dv[k].nt5 == rv[k].nt5 && dv[k].nt3 == rv[k].nt3 &&
                                         same(dv[k].cost, sv[k].cost) && same(dv[k].cost, rv[k].cost);
                                if (!ok) {
                                    if (failures < 5)
                                        std::printf("cell mismatch %s lambda=%g ML_BASE37=%d a=%d b=%d\n",
                                                    protein.c_str(), lambda, ml_base, a, b);
                                    ++failures;
                                }
                                ++cells;
                            }
                TurnerSparseStats sparse_stats, right_stats;
                sparse.sparse_stats(sparse_stats);
                right.sparse_stats(right_stats);
                if (right_stats.candidates != right_stats.direct_intervals ||
                    sparse_stats.candidates > right_stats.candidates) ++failures;
                ++cases;
            }
        }
    }
    ML_BASE37 = saved;
    DFA dfa("WIPKKSAWCFEQQLCNLYSQLFMGV", table);
    TurnerProfile sp, rp, dp;
    const double expected = turner_min_cost(dfa, 400.0);
    for (int threads : {1, 2, 4}) {
        const double s = turner_min_cost_profiled(dfa, 400.0, TurnerKernel::Sparse, sp, 3, threads);
        const double r = turner_min_cost_profiled(dfa, 400.0, TurnerKernel::DenseRightNormal, rp, 3, threads);
        const double d = turner_min_cost_profiled(dfa, 400.0, TurnerKernel::DenseWavefront, dp, 3, threads);
        if (!same(expected, s) || !same(expected, r) || !same(expected, d) ||
            sp.split_visits > rp.split_visits || sp.split_eligible > rp.split_eligible ||
            sp.multiloop_seconds < 0.0 || rp.split_visits == 0) ++failures;
    }
    std::printf("right-normal production cells: %d cases, %zu intervals, %d failures\n",
                cases, cells, failures);
    return failures ? 1 : 0;
}
