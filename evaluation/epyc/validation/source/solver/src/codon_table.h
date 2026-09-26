// codon_table.h — genetic code, codon-usage weights, and CAI.
// Clean-room: derived from the LinearDesign paper's definitions (Eqs. 3, 6) and the
// standard CAI definition (Sharp & Li 1987). No code from the restricted distribution.
#ifndef LDCLEAN_CODON_TABLE_H
#define LDCLEAN_CODON_TABLE_H

#include <map>
#include <string>
#include <vector>

namespace ldclean {

// A synonymous codon choice with its relative adaptiveness w(c) and CAI edge cost.
struct CodonOption {
    std::string codon;   // 3 nt over {A,C,G,U}
    double freq;         // raw usage frequency among synonymous codons
    double w;            // relative adaptiveness = freq / max synonymous freq, in (0,1]
    double neg_log_w;    // -log w(c)  (natural log); the per-codon CAI cost in the objective
};

// Genetic code + codon usage loaded from a CSV ("codon,aa,freq"; '*' = stop; '#' comment).
class CodonTable {
public:
    // Loads the CSV. Throws std::runtime_error on failure.
    explicit CodonTable(const std::string& csv_path);

    // Amino acids use 1-letter codes; stop is '*'.
    const std::vector<CodonOption>& options_for(char aa) const;

    // codon (ACGU) -> amino acid char ('*' for stop). 0 if unknown.
    char aa_of(const std::string& codon) const;

    // Relative adaptiveness w(c) for a codon; 0 if unknown.
    double w_of(const std::string& codon) const;

    // CAI of an mRNA coding sequence (length multiple of 3), as the geometric mean of the
    // per-codon relative adaptiveness over ALL codons in the sequence (including the stop
    // codon), matching the LinearDesign convention. Returns 0 if not divisible by 3.
    double cai(const std::string& rna) const;

    bool has_aa(char aa) const { return aa_to_options_.count(aa) > 0; }

private:
    std::map<char, std::vector<CodonOption>> aa_to_options_;   // aa -> synonymous codons
    std::map<std::string, char> codon_to_aa_;
    std::map<std::string, double> codon_to_w_;
};

}  // namespace ldclean

#endif  // LDCLEAN_CODON_TABLE_H
