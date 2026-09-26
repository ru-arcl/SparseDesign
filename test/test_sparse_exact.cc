// Exhaustive validation of the sparse weighted-lattice DP.
//
// For each small protein, enumerate every synonymous CDS explicitly. Fold each fixed CDS with the
// retained dense multiloop recurrence, add its exact codon penalty, and compare the best enumerated
// value against both the production sparse lattice solver and the dense lattice reference solver.
// This exercises three independent decompositions of the design problem:
//   (1) Cartesian-product sequence enumeration + fixed-sequence dense folding,
//   (2) one weighted codon lattice + dense folding, and
//   (3) one weighted codon lattice + exact candidate sparsification.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "../src/codon_table.h"
#include "../src/dfa.h"
#include "../src/fold_turner.h"

using namespace ldclean;

static void enumerate_cds(const CodonTable& table, const std::string& protein, size_t i,
                          std::string& sequence, std::vector<std::string>& out) {
    if (i == protein.size()) {
        out.push_back(sequence);
        return;
    }
    for (const CodonOption& option : table.options_for(protein[i])) {
        sequence += option.codon;
        enumerate_cds(table, protein, i + 1, sequence, out);
        sequence.resize(sequence.size() - 3);
    }
}

static double codon_penalty(const CodonTable& table, const std::string& rna) {
    double penalty = 0.0;
    for (size_t i = 0; i < rna.size(); i += 3)
        penalty -= std::log(table.w_of(rna.substr(i, 3)));
    return penalty;
}

int main() {
    CodonTable table("data/codon_usage_freq_table_human.csv");
    const std::vector<std::string> proteins{"MFL", "RSL", "SSSS", "MWFMWF", "MSWFWFM"};
    const std::vector<double> lambdas{0.0, 30.0, 100.0, 300.0};
    int failures = 0;
    size_t comparisons = 0;

    for (const std::string& protein : proteins) {
        std::vector<std::string> sequences;
        std::string sequence;
        enumerate_cds(table, protein, 0, sequence, sequences);

        // The fixed-sequence folding costs do not depend on lambda because the validation lattice
        // has zero edge weights. Compute each dense MFE once, then add the codon cost explicitly.
        std::vector<double> mfe(sequences.size());
        std::vector<double> penalty(sequences.size());
        for (size_t i = 0; i < sequences.size(); ++i) {
            DFA fixed(sequences[i]);
            mfe[i] = turner_min_cost_dense_reference(fixed, 0.0);
            penalty[i] = codon_penalty(table, sequences[i]);
        }

        DFA lattice(protein, table);
        for (double lambda : lambdas) {
            double brute = 1e18;
            for (size_t i = 0; i < sequences.size(); ++i)
                brute = std::min(brute, mfe[i] + lambda * penalty[i]);

            TurnerSparseStats stats;
            double sparse = turner_min_cost(lattice, lambda, 3, 1, &stats);
            double dense = turner_min_cost_dense_reference(lattice, lambda);
            bool sparse_ok = std::fabs(sparse - brute) < 1e-6;
            bool dense_ok = std::fabs(dense - brute) < 1e-6;
            if (!sparse_ok || !dense_ok) ++failures;
            ++comparisons;
            std::printf("protein=%s paths=%zu lambda=%.2f brute=%.9f dense=%.9f "
                        "sparse=%.9f Z=%zu direct=%zu %s\n",
                        protein.c_str(), sequences.size(), lambda / 100.0, brute, dense, sparse,
                        stats.candidates, stats.direct_intervals,
                        (sparse_ok && dense_ok) ? "OK" : "FAIL");
        }
    }

    // A finite infeasibility sentinel used to collide with legitimate high-lambda objectives.
    // Forbid F's maximum-weight UUC codon so the only path has a strictly positive CAI penalty;
    // its valid objective is deliberately above 1e18.
    DFA high_cost("F", table, {"UUC"});
    const double high_lambda = 1e19;
    const double high_expected = high_lambda * -std::log(table.w_of("UUU"));
    const double high_sparse = turner_min_cost(high_cost, high_lambda, 3, 4);
    const double high_dense = turner_min_cost_dense_reference(high_cost, high_lambda, 3);
    const double high_wave = turner_min_cost_dense_wavefront(high_cost, high_lambda, 3, 4);
    const double high_scale = std::fabs(high_expected);
    const bool high_ok = std::isfinite(high_sparse) && high_sparse > 1e18 &&
                         std::fabs(high_sparse - high_expected) <= 1e-12 * high_scale &&
                         std::fabs(high_dense - high_expected) <= 1e-12 * high_scale &&
                         std::fabs(high_wave - high_expected) <= 1e-12 * high_scale;
    std::printf("high-lambda expected=%.9g dense=%.9g dense_wave=%.9g sparse=%.9g %s\n",
                high_expected, high_dense, high_wave, high_sparse, high_ok ? "OK" : "FAIL");
    if (!high_ok) ++failures;
    ++comparisons;

    bool nan_rejected = false;
    bool negative_loop_rejected = false;
    try {
        (void)turner_min_cost(high_cost, std::numeric_limits<double>::quiet_NaN(), 3, 1);
    } catch (const std::invalid_argument&) {
        nan_rejected = true;
    }
    try {
        (void)turner_min_cost(high_cost, 0.0, -1, 1);
    } catch (const std::invalid_argument&) {
        negative_loop_rejected = true;
    }
    const bool argument_ok = nan_rejected && negative_loop_rejected;
    std::printf("invalid-argument nan=%s negative_min_loop=%s %s\n",
                nan_rejected ? "rejected" : "accepted",
                negative_loop_rejected ? "rejected" : "accepted",
                argument_ok ? "OK" : "FAIL");
    if (!argument_ok) ++failures;
    ++comparisons;

    // A representational overflow is a valid +infinity result for the cost-only API, but the
    // design API must stop before decision replay rather than backtracing an absent finite path.
    DFA overflow_cost("L", table, {"UUG", "CUU", "CUC", "CUA", "CUG"});
    const double overflow_lambda = std::numeric_limits<double>::max();
    const bool cost_overflowed = std::isinf(turner_min_cost(overflow_cost, overflow_lambda, 3, 1));
    bool design_overflow_rejected = false;
    try {
        (void)design_turner(overflow_cost, table, overflow_lambda, 3, 1);
    } catch (const std::overflow_error&) {
        design_overflow_rejected = true;
    }
    const bool overflow_ok = cost_overflowed && design_overflow_rejected;
    std::printf("overflow cost=%s design=%s %s\n",
                cost_overflowed ? "+inf" : "finite",
                design_overflow_rejected ? "rejected" : "backtraced",
                overflow_ok ? "OK" : "FAIL");
    if (!overflow_ok) ++failures;
    ++comparisons;

    std::printf("exhaustive sparse checks: %zu comparisons, %d failures\n", comparisons, failures);
    return failures == 0 ? 0 : 1;
}
