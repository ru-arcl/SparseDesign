// test_dfa.cc — verify the weighted DFA encodes exactly the synonymous design space,
// with path weights == sum of -log w(codon) (so exp(-weight/#codons) == CAI).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <set>
#include <string>
#include <vector>

#include "../src/codon_table.h"
#include "../src/dfa.h"

using namespace ldclean;

// Enumerate all start->final paths: collect (sequence, total edge weight).
static void enumerate(const DFA& d, int node, std::string& seq, double w,
                      std::vector<std::pair<std::string, double>>& out) {
    if (node == d.final()) { out.emplace_back(seq, w); return; }
    for (const Edge& e : d.out_edges(node)) {
        seq.push_back(nuc_char(e.nuc));
        enumerate(d, e.to, seq, w + e.weight, out);
        seq.pop_back();
    }
}

// Brute-force the synonymous design space (Cartesian product of codon choices).
static void product(const CodonTable& t, const std::string& residues, size_t k,
                    std::string& acc, std::set<std::string>& out) {
    if (k == residues.size()) { out.insert(acc); return; }
    for (const auto& opt : t.options_for(residues[k])) {
        acc += opt.codon;
        product(t, residues, k + 1, acc, out);
        acc.erase(acc.size() - opt.codon.size());
    }
}

// Position-major node IDs are part of the compact Turner memo-layout contract. Also verify that
// relabeling kept the layered edge endpoints internally consistent in both adjacency directions.
static bool check_position_major(const DFA& d, const char* label) {
    int expected_node = 0;
    int bad_nodes = 0;
    int bad_edges = 0;
    for (int pos = 0; pos <= d.rna_length(); ++pos)
        for (int node : d.nodes_at(pos)) {
            if (node != expected_node || d.pos_of(node) != pos) ++bad_nodes;
            ++expected_node;
        }
    if (expected_node != d.num_nodes()) ++bad_nodes;

    for (int node = 0; node < d.num_nodes(); ++node) {
        for (const Edge& edge : d.out_edges(node))
            if (edge.from != node || edge.to < 0 || edge.to >= d.num_nodes() ||
                d.pos_of(edge.to) != d.pos_of(node) + 1)
                ++bad_edges;
        for (const Edge& edge : d.in_edges(node))
            if (edge.to != node || edge.from < 0 || edge.from >= d.num_nodes() ||
                d.pos_of(edge.from) + 1 != d.pos_of(node))
                ++bad_edges;
    }

    const bool endpoints_ok = d.start() == 0 && d.final() == d.num_nodes() - 1;
    const bool ok = bad_nodes == 0 && bad_edges == 0 && endpoints_ok;
    std::printf("position-major %-11s : %s (bad_nodes=%d bad_edges=%d endpoints=%s)\n",
                label, ok ? "OK" : "FAIL", bad_nodes, bad_edges,
                endpoints_ok ? "OK" : "FAIL");
    return ok;
}

int main() {
    CodonTable t("data/codon_usage_freq_table_human.csv");
    int failures = 0;

    // "MFL" exercises a 1-codon (M), 2-codon (F), and 6-codon (L) amino acid.
    std::string protein = "MFL";
    DFA d(protein, t);
    if (!check_position_major(d, "synonymous")) ++failures;

    std::vector<std::pair<std::string, double>> paths;
    std::string seq;
    enumerate(d, d.start(), seq, 0.0, paths);

    std::set<std::string> dfa_seqs;
    int weight_mismatch = 0;
    for (auto& pw : paths) {
        dfa_seqs.insert(pw.first);
        int ncodons = static_cast<int>(pw.first.size() / 3);
        double cai_from_weight = std::exp(-pw.second / ncodons);
        if (std::fabs(cai_from_weight - t.cai(pw.first)) > 1e-9) ++weight_mismatch;
    }

    std::set<std::string> expected;
    std::string acc;
    std::string residues = protein;  // LinearDesign designs exactly the input residues (no auto-stop)
    product(t, residues, 0, acc, expected);

    // Expected count = product of synonymous-codon counts (M=1, F=2, L=6) = 12.
    printf("DFA paths: %zu   distinct DFA seqs: %zu   expected seqs: %zu\n",
           paths.size(), dfa_seqs.size(), expected.size());
    printf("nodes: %d   rna_length: %d (expected %d)\n", d.num_nodes(), d.rna_length(),
           3 * (int)residues.size());

    bool langs_equal = (dfa_seqs == expected);
    printf("language == synonymous space : %s\n", langs_equal ? "OK" : "FAIL");
    if (!langs_equal) ++failures;
    printf("no duplicate paths           : %s\n",
           paths.size() == dfa_seqs.size() ? "OK" : "FAIL");
    if (paths.size() != dfa_seqs.size()) ++failures;
    printf("path weight == -L*log CAI    : %s (%d mismatches)\n",
           weight_mismatch == 0 ? "OK" : "FAIL", weight_mismatch);
    if (weight_mismatch != 0) ++failures;

    // Every sequence must translate back to the protein.
    int bad_translation = 0;
    for (const auto& s : dfa_seqs) {
        std::string aas;
        for (size_t i = 0; i < s.size(); i += 3) aas.push_back(t.aa_of(s.substr(i, 3)));
        if (aas != residues) ++bad_translation;
    }
    printf("all paths translate to MFL   : %s (%d bad)\n",
           bad_translation == 0 ? "OK" : "FAIL", bad_translation);
    if (bad_translation != 0) ++failures;

    // Larger protein: just check counts line up with the product of option counts.
    std::string big = "MYGKIIFVLLLSGIVS";
    DFA d2(big, t);
    if (!check_position_major(d2, "large")) ++failures;
    double expected_paths = 1.0;
    for (char aa : big) expected_paths *= t.options_for(aa).size();
    printf("\nbig protein nodes=%d rna_len=%d ; #designs=%.0f (<= 2n nodes? %d<=%d %s)\n",
           d2.num_nodes(), d2.rna_length(), expected_paths, d2.num_nodes(),
           2 * d2.rna_length(), d2.num_nodes() <= 2 * d2.rna_length() ? "OK" : "FAIL");
    if (d2.num_nodes() > 2 * d2.rna_length()) ++failures;

    // Cover the other public construction paths, including multiple product states on one layer
    // and the zero-length endpoint alias.
    DFA fixed(std::string("AUGU"));
    DFA constrained("SSSS", t, {"UCUUCU", "AGCAGC"});
    DFA empty(std::string(""));
    if (!check_position_major(fixed, "fixed")) ++failures;
    if (!check_position_major(constrained, "product")) ++failures;
    if (!check_position_major(empty, "empty")) ++failures;

    printf("\n%s\n", failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED");
    return failures == 0 ? 0 : 1;
}
