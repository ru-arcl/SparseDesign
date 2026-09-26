// dfa.cc — see dfa.h.
#include "dfa.h"

#include <algorithm>
#include <cctype>
#include <map>
#include <stdexcept>
#include <utility>

namespace ldclean {

int DFA::add_node(int pos) {
    int id = static_cast<int>(out_.size());
    out_.emplace_back();
    in_.emplace_back();
    node_pos_.push_back(pos);
    if (pos >= static_cast<int>(by_pos_.size())) by_pos_.resize(pos + 1);
    by_pos_[pos].push_back(id);
    return id;
}

void DFA::add_edge(int from, int nuc, double weight, int to) {
    // De-duplicate identical edges (shared codon prefixes create the same prefix edge).
    for (const auto& e : out_[from])
        if (e.nuc == nuc && e.to == to) return;
    out_[from].push_back(Edge{nuc, weight, to, from});
    in_[to].push_back(Edge{nuc, weight, to, from});
}

// Stable post-build relabeling: node IDs become the concatenation of nodes_at(0), nodes_at(1),
// ... without changing any within-layer or adjacency-vector order. The folding DP can then use a
// node ID directly as its position-major rank while preserving every recurrence's tie order.
void DFA::renumber_position_major() {
    const int n = num_nodes();
    std::vector<int> old_to_new(n, -1);
    int next = 0;
    bool identity = true;
    for (int pos = 0; pos < static_cast<int>(by_pos_.size()); ++pos)
        for (int old : by_pos_[pos]) {
            if (old < 0 || old >= n || old_to_new[old] != -1 || node_pos_[old] != pos)
                throw std::logic_error("DFA layers do not partition node IDs by position");
            old_to_new[old] = next;
            identity = identity && old == next;
            ++next;
        }
    if (next != n) throw std::logic_error("DFA layers do not cover every node ID");
    if (identity) return;

    std::vector<int> new_node_pos(n);
    std::vector<std::vector<Edge>> new_out(n), new_in(n);
    std::vector<std::vector<int>> new_by_pos(by_pos_.size());
    for (int pos = 0; pos < static_cast<int>(by_pos_.size()); ++pos) {
        new_by_pos[pos].reserve(by_pos_[pos].size());
        for (int old : by_pos_[pos]) {
            const int now = old_to_new[old];
            new_node_pos[now] = pos;
            new_out[now] = std::move(out_[old]);
            new_in[now] = std::move(in_[old]);
            new_by_pos[pos].push_back(now);
        }
    }
    for (int node = 0; node < n; ++node) {
        for (Edge& edge : new_out[node]) {
            edge.from = old_to_new[edge.from];
            edge.to = old_to_new[edge.to];
        }
        for (Edge& edge : new_in[node]) {
            edge.from = old_to_new[edge.from];
            edge.to = old_to_new[edge.to];
        }
    }
    start_ = old_to_new[start_];
    final_ = old_to_new[final_];
    node_pos_ = std::move(new_node_pos);
    out_ = std::move(new_out);
    in_ = std::move(new_in);
    by_pos_ = std::move(new_by_pos);
}

DFA::DFA(const std::string& protein, const CodonTable& table) {
    // Design exactly the input residues (LinearDesign does NOT auto-append a stop codon;
    // include '*' in the protein if a stop codon is wanted).
    const std::string& residues = protein;

    start_ = add_node(0);
    int boundary = start_;  // start node of the current codon (always at a codon boundary)

    for (char aa : residues) {
        if (!table.has_aa(aa))
            throw std::runtime_error(std::string("protein contains unknown amino acid: ") + aa);
        const auto& options = table.options_for(aa);
        int base_pos = node_pos_[boundary];
        int end_node = add_node(base_pos + 3);  // shared codon-boundary end node

        // Trie keyed by codon prefix ("" = boundary); intermediate nodes shared across codons.
        std::map<std::string, int> prefix_node;
        prefix_node[""] = boundary;

        for (const auto& opt : options) {
            const std::string& codon = opt.codon;  // 3 nt
            int cur = boundary;
            for (int d = 0; d < 3; ++d) {
                int nuc = nuc_code(codon[d]);
                if (nuc == 0) throw std::runtime_error("bad nucleotide in codon " + codon);
                if (d < 2) {
                    std::string key = codon.substr(0, d + 1);
                    auto it = prefix_node.find(key);
                    int nxt;
                    if (it == prefix_node.end()) {
                        nxt = add_node(base_pos + d + 1);
                        prefix_node[key] = nxt;
                    } else {
                        nxt = it->second;
                    }
                    add_edge(cur, nuc, 0.0, nxt);
                    cur = nxt;
                } else {
                    // third nucleotide carries this codon's CAI cost; converges to end_node
                    add_edge(cur, nuc, opt.neg_log_w, end_node);
                }
            }
        }
        boundary = end_node;
    }

    final_ = boundary;
    rna_length_ = node_pos_[final_];
    renumber_position_major();
}

DFA::DFA(const std::string& rna) {
    start_ = add_node(0);
    int cur = start_;
    for (size_t i = 0; i < rna.size(); ++i) {
        int nuc = nuc_code(rna[i]);
        if (nuc == 0)
            throw std::runtime_error(std::string("fixed RNA contains bad nucleotide: ") + rna[i]);
        int next = add_node(static_cast<int>(i) + 1);
        add_edge(cur, nuc, 0.0, next);
        cur = next;
    }
    final_ = cur;
    rna_length_ = static_cast<int>(rna.size());
    renumber_position_major();
}

namespace {

bool ends_with(const std::string& text, const std::string& suffix) {
    return suffix.size() <= text.size() &&
           text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// Small deterministic matcher represented by all motif-prefix suffix states. This is equivalent
// to an Aho--Corasick automaton; constructing its four transitions by direct suffix tests keeps the
// validation-facing implementation compact and makes the accepted language easy to audit.
class ForbiddenMotifDFA {
public:
    explicit ForbiddenMotifDFA(std::vector<std::string> motifs) {
        for (std::string& motif : motifs) {
            if (motif.empty()) throw std::runtime_error("forbidden motif must not be empty");
            for (char& ch : motif) {
                ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
                if (ch == 'T') ch = 'U';
                if (nuc_code(ch) == 0)
                    throw std::runtime_error("forbidden motifs must contain only A, C, G, U/T");
            }
        }
        std::sort(motifs.begin(), motifs.end());
        motifs.erase(std::unique(motifs.begin(), motifs.end()), motifs.end());
        motifs_ = std::move(motifs);

        prefixes_.push_back("");
        for (const std::string& motif : motifs_)
            for (size_t n = 1; n < motif.size(); ++n)
                prefixes_.push_back(motif.substr(0, n));
        std::sort(prefixes_.begin(), prefixes_.end(), [](const std::string& a,
                                                         const std::string& b) {
            return a.size() != b.size() ? a.size() < b.size() : a < b;
        });
        prefixes_.erase(std::unique(prefixes_.begin(), prefixes_.end()), prefixes_.end());

        next_.assign(prefixes_.size(), std::vector<int>(5, -1));
        for (size_t state = 0; state < prefixes_.size(); ++state)
            for (int nuc = 1; nuc <= 4; ++nuc) {
                std::string extended = prefixes_[state] + nuc_char(nuc);
                bool forbidden = false;
                for (const std::string& motif : motifs_)
                    if (ends_with(extended, motif)) { forbidden = true; break; }
                if (forbidden) continue;

                int best = 0;
                for (size_t candidate = 1; candidate < prefixes_.size(); ++candidate)
                    if (prefixes_[candidate].size() > prefixes_[best].size() &&
                        ends_with(extended, prefixes_[candidate]))
                        best = static_cast<int>(candidate);
                next_[state][nuc] = best;
            }
    }

    int next(int state, int nuc) const { return next_[state][nuc]; }

private:
    std::vector<std::string> motifs_;
    std::vector<std::string> prefixes_;
    std::vector<std::vector<int>> next_;
};

}  // namespace

DFA::DFA(const std::string& protein, const CodonTable& table,
         const std::vector<std::string>& forbidden_motifs) {
    DFA base(protein, table);
    if (forbidden_motifs.empty()) {
        node_pos_ = std::move(base.node_pos_);
        out_ = std::move(base.out_);
        in_ = std::move(base.in_);
        by_pos_ = std::move(base.by_pos_);
        start_ = base.start_;
        final_ = base.final_;
        rna_length_ = base.rna_length_;
        return;
    }

    ForbiddenMotifDFA constraint(forbidden_motifs);
    rna_length_ = base.rna_length_;
    start_ = add_node(0);
    if (rna_length_ == 0) { final_ = start_; return; }
    // All accepting matcher states can merge at the sequence endpoint because no future symbol
    // can distinguish them. Intermediate states retain matcher context across codon boundaries.
    final_ = add_node(rna_length_);
    std::vector<std::map<int, int>> product(base.num_nodes());
    product[base.start()][0] = start_;

    for (int pos = 0; pos < rna_length_; ++pos)
        for (int base_from : base.nodes_at(pos))
            for (const auto& state_node : product[base_from]) {
                int matcher_state = state_node.first;
                int from = state_node.second;
                for (const Edge& edge : base.out_edges(base_from)) {
                    int next_state = constraint.next(matcher_state, edge.nuc);
                    if (next_state < 0) continue;
                    int to;
                    if (edge.to == base.final()) {
                        to = final_;
                    } else {
                        auto found = product[edge.to].find(next_state);
                        if (found == product[edge.to].end()) {
                            to = add_node(base.pos_of(edge.to));
                            product[edge.to][next_state] = to;
                        } else {
                            to = found->second;
                        }
                    }
                    add_edge(from, edge.nuc, edge.weight, to);
                }
            }

    if (in_[final_].empty())
        throw std::runtime_error("forbidden motifs eliminate every synonymous coding sequence");
    renumber_position_major();
}

}  // namespace ldclean
