// dfa.h — weighted codon DFA over the mRNA design space.
// Clean-room: derived from the LinearDesign paper ("DFA representations for codons and mRNA
// candidate sequences"; Sup Figs 5-7). Nucleotides label edges; the per-codon CAI cost
// -log w(c) is placed on each codon's third-nucleotide edge so every start->final path's
// total edge weight equals -log CAI(path) * (#codons) = sum of -log w over its codons.
#ifndef LDCLEAN_DFA_H
#define LDCLEAN_DFA_H

#include <string>
#include <vector>

#include "codon_table.h"

namespace ldclean {

// Nucleotide encoding (Vienna): A=1, C=2, G=3, U=4.
inline int nuc_code(char c) {
    switch (c) { case 'A': return 1; case 'C': return 2; case 'G': return 3; case 'U': return 4; }
    return 0;
}
inline char nuc_char(int code) {
    switch (code) { case 1: return 'A'; case 2: return 'C'; case 3: return 'G'; case 4: return 'U'; }
    return 'N';
}

struct Edge {
    int nuc;        // 1..4 (A,C,G,U)
    double weight;  // CAI edge cost (-log w on a codon's 3rd-nt edge; else 0)
    int to;         // target node id
    int from;       // source node id (used for in-edges)
};

// A weighted DFA whose start->final paths are exactly the synonymous coding sequences of a target
// protein, or a single fixed RNA path when constructed through the validation overload.
class DFA {
public:
    // Build D(p0) o D(p1) o ... for exactly the supplied residues. Include '*' explicitly when
    // a stop codon is desired.
    DFA(const std::string& protein, const CodonTable& table);

    // Build a single-path, zero-weight lattice spelling exactly `rna`. This is primarily a
    // validation hook: exhaustive synonymous-sequence tests fold each fixed sequence with the
    // retained dense reference recurrence, independently of the sparse lattice recurrence.
    explicit DFA(const std::string& rna);

    // Compose the synonymous-codon lattice with a deterministic finite-state constraint that
    // rejects every path containing any forbidden RNA motif. Motifs may cross codon boundaries.
    // The product remains layered and acyclic, so the folding solver is unchanged.
    DFA(const std::string& protein, const CodonTable& table,
        const std::vector<std::string>& forbidden_motifs);

    int num_nodes() const { return static_cast<int>(out_.size()); }
    int start() const { return start_; }
    int final() const { return final_; }
    int rna_length() const { return rna_length_; }            // = 3 * (#codons incl. stop)
    int pos_of(int node) const { return node_pos_[node]; }

    const std::vector<Edge>& out_edges(int node) const { return out_[node]; }
    const std::vector<Edge>& in_edges(int node) const { return in_[node]; }

    // Node ids grouped by nucleotide position (0..rna_length). nodes_at(i) = ids at position i.
    const std::vector<int>& nodes_at(int pos) const { return by_pos_[pos]; }

private:
    int add_node(int pos);
    void add_edge(int from, int nuc, double weight, int to);
    void renumber_position_major();

    std::vector<int> node_pos_;
    std::vector<std::vector<Edge>> out_;
    std::vector<std::vector<Edge>> in_;
    std::vector<std::vector<int>> by_pos_;
    int start_ = 0;
    int final_ = 0;
    int rna_length_ = 0;
};

}  // namespace ldclean

#endif  // LDCLEAN_DFA_H
