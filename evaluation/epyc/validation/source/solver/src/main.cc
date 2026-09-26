// main.cc — clean-room LinearDesign CLI.
// Reads a protein (FASTA or raw, from stdin), designs an optimal mRNA CDS jointly minimizing
// MFE and maximizing CAI, and prints the CDS, its MFE secondary structure, the MFE, and CAI.
//
//   echo MFLMVF | ./lineardesign-clean [-l LAMBDA] [-c CODON_TABLE]
//   ./lineardesign-clean -l 1.0 < protein.fasta
//   echo MFLMVF | ./lineardesign-clean --evaluate-mrna AUGUUCCUGAUGGUGUUC
//
// LAMBDA is in the LinearDesign paper's units (0 = pure MFE; larger favors CAI); internally
// it is scaled by 100 to combine with integer (0.01 kcal/mol) energies. Designs exactly the
// input residues (no auto-appended stop codon; include '*' in the protein for a stop).
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "codon_table.h"
#include "dfa.h"
#include "fold_turner.h"

namespace {

void print_usage(const char* program) {
    fprintf(stderr, "usage: %s [-l LAMBDA] [-c CODON_TABLE] [-j THREADS] [-v] "
                    "[--forbid-motif RNA] [--sparse-stats] [--diagnostics] [--evaluate-mrna RNA] "
                    " (protein on stdin)\n", program);
}

std::string normalize_rna(const std::string& text) {
    std::string rna;
    for (char raw : text) {
        unsigned char byte = static_cast<unsigned char>(raw);
        if (std::isspace(byte)) continue;
        char ch = static_cast<char>(std::toupper(byte));
        if (ch == 'T') ch = 'U';
        if (ldclean::nuc_code(ch) == 0)
            throw std::invalid_argument("--evaluate-mrna expects only A, C, G, U/T");
        rna.push_back(ch);
    }
    if (rna.empty()) throw std::invalid_argument("--evaluate-mrna requires a nonempty RNA");
    return rna;
}

double parse_lambda(const std::string& text) {
    size_t consumed = 0;
    double value = 0.0;
    try {
        value = std::stod(text, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument("--lambda expects a finite nonnegative number: '" + text + "'");
    }
    if (consumed != text.size() || !std::isfinite(value) || value < 0.0 ||
        value > std::numeric_limits<double>::max() / 100.0)
        throw std::invalid_argument("--lambda expects a finite nonnegative number: '" + text + "'");
    return value;
}

int parse_threads(const std::string& text) {
    size_t consumed = 0;
    long long value = 0;
    try {
        value = std::stoll(text, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument("--threads expects a positive integer: '" + text + "'");
    }
    if (consumed != text.size() || value < 1 || value > std::numeric_limits<int>::max())
        throw std::invalid_argument("--threads expects a positive integer: '" + text + "'");
    return static_cast<int>(value);
}

int run(int argc, char** argv) {
    double lambda_user = 0.0;
    std::string codon_table = "data/codon_usage_freq_table_human.csv";
    std::string evaluate_mrna;
    std::vector<std::string> forbidden_motifs;
    bool verbose = false;
    bool sparse_stats = false;
    bool diagnostics = false;
    int threads = 1;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-l" || a == "--lambda") {
            if (i + 1 >= argc) throw std::invalid_argument("missing value for " + a);
            lambda_user = parse_lambda(argv[++i]);
        } else if (a == "-c" || a == "--codonusage") {
            if (i + 1 >= argc) throw std::invalid_argument("missing value for " + a);
            codon_table = argv[++i];
            if (codon_table.empty()) throw std::invalid_argument(a + " requires a nonempty path");
        } else if (a == "-j" || a == "--threads") {
            if (i + 1 >= argc) throw std::invalid_argument("missing value for " + a);
            threads = parse_threads(argv[++i]);
        } else if (a == "--forbid-motif") {
            if (i + 1 >= argc) throw std::invalid_argument("missing value for " + a);
            forbidden_motifs.push_back(argv[++i]);
            if (forbidden_motifs.back().empty())
                throw std::invalid_argument("--forbid-motif requires a nonempty RNA motif");
        } else if (a == "--evaluate-mrna") {
            if (i + 1 >= argc) throw std::invalid_argument("missing value for " + a);
            evaluate_mrna = normalize_rna(argv[++i]);
        }
        else if (a == "-v" || a == "--verbose") verbose = true;
        else if (a == "--sparse-stats") sparse_stats = true;
        else if (a == "--diagnostics") diagnostics = true;
        else if (a == "-h" || a == "--help") {
            print_usage(argv[0]);
            return 0;
        } else throw std::invalid_argument("unknown option: " + a);
    }

    // One protein per invocation. A second FASTA header must never concatenate unrelated genes.
    std::string name = ">designed", protein;
    bool seen_header = false;
    std::string line;
    while (std::getline(std::cin, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n' || isspace((unsigned char)line.back())))
            line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (seen_header || !protein.empty())
                throw std::invalid_argument("expected one protein record; multiple FASTA records or mixed raw/FASTA input");
            seen_header = true;
            name = line;
            continue;
        }
        for (char ch : line) if (!isspace((unsigned char)ch)) protein.push_back(toupper((unsigned char)ch));
    }
    if (protein.empty()) { fprintf(stderr, "error: no protein on stdin\n"); return 1; }

    ldclean::CodonTable table(codon_table);
    for (char aa : protein)
        if (!table.has_aa(aa)) { fprintf(stderr, "error: unknown amino acid '%c'\n", aa); return 1; }

    ldclean::TurnerSparseStats stats;
    ldclean::TurnerDesign r;
    if (!evaluate_mrna.empty()) {
        if (lambda_user != 0.0 || !forbidden_motifs.empty() || sparse_stats)
            throw std::invalid_argument(
                "--evaluate-mrna cannot be combined with nonzero --lambda, --forbid-motif, "
                "or --sparse-stats");
        if (evaluate_mrna.size() != protein.size() * 3)
            throw std::invalid_argument("evaluated mRNA length does not equal three times the protein length");
        for (size_t i = 0; i < protein.size(); ++i) {
            std::string codon = evaluate_mrna.substr(3 * i, 3);
            if (table.aa_of(codon) != protein[i])
                throw std::invalid_argument("evaluated mRNA does not encode the input protein at codon " +
                                            std::to_string(i + 1) + " (" + codon + ")");
        }
        // A one-path, zero-weight DFA computes the supplied sequence's MFE structure.  Its CAI is
        // evaluated independently from the selected codons below.
        ldclean::DFA fixed(evaluate_mrna);
        r = ldclean::design_turner(fixed, table, 0.0, 3, threads);
    } else {
        std::unique_ptr<ldclean::DFA> dfa;
        if (forbidden_motifs.empty()) dfa.reset(new ldclean::DFA(protein, table));
        else dfa.reset(new ldclean::DFA(protein, table, forbidden_motifs));
        r = ldclean::design_turner(
            *dfa, table, lambda_user * 100.0, 3, threads, sparse_stats ? &stats : nullptr);
    }

    std::cout << name << "\n";
    if (verbose) { std::cout << "Input protein: " << protein << "\n";
                   int eff = ldclean::turner_effective_threads(threads);
                   std::cout << "Threads: " << eff;
                   if (eff != threads) std::cout << " (requested " << threads
                                                 << ", capped to host cores)";
                   std::cout << "\n"; }
    std::cout << "mRNA sequence:  " << r.mrna << "\n";
    std::cout << "mRNA structure: " << r.structure << "\n";
    printf("mRNA folding free energy: %.2f kcal/mol; mRNA CAI: %.3f\n", r.mfe_kcal, r.cai);
    if (diagnostics)
        fprintf(stderr, "diagnostics dp_cost=%.17g realized_cost=%.17g objective_delta=%.17g "
                        "codon_penalty=%.17g cai=%.17g cost_unit=centikcal_mol\n",
                r.cost, r.realized_cost, r.realized_cost - r.cost, r.codon_penalty, r.cai);
    if (sparse_stats) {
        double direct_fraction = stats.direct_intervals
            ? static_cast<double>(stats.candidates) / stats.direct_intervals : 0.0;
        double quadratic_fraction = stats.rna_length
            ? static_cast<double>(stats.candidates) /
              (static_cast<double>(stats.rna_length) * stats.rna_length) : 0.0;
        fprintf(stderr,
                "sparse_stats rna_nt=%d lattice_nodes=%zu intervals=%zu direct=%zu "
                "candidates=%zu candidate_over_direct=%.9g Z_over_n2=%.9g "
                "capacity=%zu storage_bytes=%zu max_per_end=%zu\n",
                stats.rna_length, stats.lattice_nodes, stats.interval_cells,
                stats.direct_intervals, stats.candidates, direct_fraction,
                quadratic_fraction, stats.candidate_capacity, stats.candidate_storage_bytes,
                stats.max_candidates_per_end);
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run(argc, argv);
    } catch (const std::exception& error) {
        fprintf(stderr, "error: %s\n", error.what());
        return 1;
    } catch (...) {
        fprintf(stderr, "error: unexpected failure\n");
        return 1;
    }
}
