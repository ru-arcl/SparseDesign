// Validate exact composition of the synonymous-codon lattice with forbidden-motif DFAs.
// Motifs include a DNA-alphabet spelling (TAG -> UAG) and cases that occur across codon
// boundaries, so filtering codons independently would not pass this test.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "../src/codon_table.h"
#include "../src/dfa.h"
#include "../src/fold_turner.h"

using namespace ldclean;

// Test-only access to the vendored Turner parameter table. The production size-zero hairpin entry
// is prohibitive, which can mask an invalid empty DFA gap behind an energy that never wins. Making
// that one legal test parameter attractive turns the language error into an observable optimum;
// it is restored before exit.
extern int hairpin37[31];

static void enumerate_cds(const CodonTable& table, const std::string& protein, size_t i,
                          std::string& sequence, std::vector<std::string>& out) {
    if (i == protein.size()) { out.push_back(sequence); return; }
    for (const CodonOption& option : table.options_for(protein[i])) {
        sequence += option.codon;
        enumerate_cds(table, protein, i + 1, sequence, out);
        sequence.resize(sequence.size() - 3);
    }
}

static void enumerate_paths(const DFA& dfa, int node, std::string& sequence, double weight,
                            std::map<std::string, double>& out, int& duplicate_paths) {
    if (node == dfa.final()) {
        if (!out.emplace(sequence, weight).second) ++duplicate_paths;
        return;
    }
    for (const Edge& edge : dfa.out_edges(node)) {
        sequence.push_back(nuc_char(edge.nuc));
        enumerate_paths(dfa, edge.to, sequence, weight + edge.weight, out, duplicate_paths);
        sequence.pop_back();
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
    const std::string protein = "SSSS";
    const std::vector<std::string> input_motifs{"TAG", "UCG", "CAG"};
    const std::vector<std::string> rna_motifs{"UAG", "UCG", "CAG"};
    int failures = 0;

    std::vector<std::string> all;
    std::string sequence;
    enumerate_cds(table, protein, 0, sequence, all);
    std::set<std::string> expected;
    std::vector<int> removed_by(rna_motifs.size(), 0);
    for (const std::string& rna : all) {
        bool forbidden = false;
        for (size_t i = 0; i < rna_motifs.size(); ++i)
            if (rna.find(rna_motifs[i]) != std::string::npos) {
                forbidden = true;
                ++removed_by[i];
            }
        if (!forbidden) expected.insert(rna);
    }

    DFA constrained(protein, table, input_motifs);
    std::map<std::string, double> actual;
    int duplicate_paths = 0;
    enumerate_paths(constrained, constrained.start(), sequence, 0.0, actual, duplicate_paths);
    std::set<std::string> actual_language;
    int bad_weight = 0;
    for (const auto& path : actual) {
        actual_language.insert(path.first);
        if (std::fabs(path.second - codon_penalty(table, path.first)) > 1e-9) ++bad_weight;
    }
    bool every_motif_exercised =
        std::find(removed_by.begin(), removed_by.end(), 0) == removed_by.end();
    bool language_ok = actual_language == expected;
    if (!language_ok || duplicate_paths || bad_weight || !every_motif_exercised) ++failures;
    std::printf("motif product: base=%zu expected=%zu actual=%zu nodes=%d duplicates=%d "
                "bad_weights=%d each_motif_exercised=%s %s\n",
                all.size(), expected.size(), actual.size(), constrained.num_nodes(), duplicate_paths,
                bad_weight, every_motif_exercised ? "yes" : "no", language_ok ? "OK" : "FAIL");

    std::vector<double> mfe;
    std::vector<double> penalty;
    mfe.reserve(expected.size());
    penalty.reserve(expected.size());
    for (const std::string& rna : expected) {
        DFA fixed(rna);
        mfe.push_back(turner_min_cost_dense_reference(fixed, 0.0));
        penalty.push_back(codon_penalty(table, rna));
    }

    const std::vector<double> lambdas{0.0, 30.0, 100.0, 300.0};
    for (double lambda : lambdas) {
        double brute = 1e18;
        for (size_t i = 0; i < mfe.size(); ++i)
            brute = std::min(brute, mfe[i] + lambda * penalty[i]);
        TurnerSparseStats stats;
        double sparse = turner_min_cost(constrained, lambda, 3, 1, &stats);
        double dense = turner_min_cost_dense_reference(constrained, lambda);
        bool ok = std::fabs(sparse - brute) < 1e-6 && std::fabs(dense - brute) < 1e-6;
        if (!ok) ++failures;
        std::printf("motif optimum lambda=%.2f brute=%.9f dense=%.9f sparse=%.9f Z=%zu %s\n",
                    lambda / 100.0, brute, dense, sparse, stats.candidates, ok ? "OK" : "FAIL");
    }

    // With min_loop=0, a size-zero hairpin asks gap_options() about the two inner endpoints on
    // the same nucleotide layer. Distinct product-DFA states are not epsilon-connected, so an
    // empty gap is legal only when those endpoints are the same node. Compare against explicit
    // fixed-sequence folding, whose one-state-per-layer DFA provides an independent oracle for
    // this otherwise easy-to-hide product-lattice case.
    const std::vector<std::string> zero_loop_motifs{"AA"};
    std::vector<std::string> zero_loop_language;
    for (const std::string& rna : all)
        if (rna.find("AA") == std::string::npos) zero_loop_language.push_back(rna);
    DFA zero_loop_constrained(protein, table, zero_loop_motifs);
    const int saved_hairpin0 = hairpin37[0];
    hairpin37[0] = -1000;
    const std::vector<double> zero_loop_lambdas{0.0, 100.0};
    for (double lambda : zero_loop_lambdas) {
        double brute = 1e18;
        for (const std::string& rna : zero_loop_language) {
            DFA fixed(rna);
            const double value = turner_min_cost_dense_reference(fixed, 0.0, 0) +
                                 lambda * codon_penalty(table, rna);
            brute = std::min(brute, value);
        }
        const double sparse = turner_min_cost(zero_loop_constrained, lambda, 0, 1);
        const double dense = turner_min_cost_dense_reference(zero_loop_constrained, lambda, 0);
        const double wave = turner_min_cost_dense_wavefront(zero_loop_constrained, lambda, 0, 4);
        const bool ok = std::fabs(sparse - brute) < 1e-6 &&
                        std::fabs(dense - brute) < 1e-6 &&
                        std::fabs(wave - brute) < 1e-6;
        if (!ok) ++failures;
        std::printf("zero-loop product lambda=%.2f paths=%zu brute=%.9f dense=%.9f "
                    "dense_wave_j4=%.9f sparse=%.9f %s\n",
                    lambda / 100.0, zero_loop_language.size(), brute, dense, wave, sparse,
                    ok ? "OK" : "FAIL");
    }
    hairpin37[0] = saved_hairpin0;

    std::printf("forbidden-motif checks: %d failures\n", failures);
    return failures == 0 ? 0 : 1;
}
