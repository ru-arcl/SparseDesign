// Three-arm, cost-only publication harness; see publication_campaign.py.
#include <chrono>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include "codon_table.h"
#include "dfa.h"
#include "fold_turner.h"

int main(int argc, char** argv) {
    try {
        if (argc < 6 || argc > 7)
            throw std::invalid_argument("usage: publication-bench MODE LAMBDA THREADS CODON_TABLE FASTA [profile]");
        const std::string mode(argv[1]);
        const double lambda = std::stod(argv[2]);
        const int threads = std::stoi(argv[3]);
        if (!std::isfinite(lambda) || lambda < 0 || threads < 1)
            throw std::invalid_argument("invalid lambda or threads");
        ldclean::TurnerKernel kernel;
        if (mode == "sparse") kernel = ldclean::TurnerKernel::Sparse;
        else if (mode == "dense_scalar") kernel = ldclean::TurnerKernel::DenseRightNormal;
        else if (mode == "dense_typed") kernel = ldclean::TurnerKernel::DenseWavefront;
        else throw std::invalid_argument("unknown kernel");
        const bool profile = argc == 7 && std::string(argv[6]) == "profile";
        if (argc == 7 && !profile) throw std::invalid_argument("unknown final argument");
        std::ifstream input(argv[5]);
        if (!input) throw std::invalid_argument("cannot open FASTA");
        std::string sequence, line;
        int headers = 0;
        while (std::getline(input, line)) {
            if (!line.empty() && line[0] == '>') {
                if (++headers > 1) throw std::invalid_argument("multiple FASTA records");
            } else for (char c : line) {
                if (!std::isspace(static_cast<unsigned char>(c))) sequence += std::toupper(static_cast<unsigned char>(c));
            }
        }
        if (sequence.empty()) throw std::invalid_argument("empty protein");
        ldclean::CodonTable table(argv[4]);
        ldclean::DFA dfa(sequence, table);
        ldclean::TurnerProfile counters;
        const int effective = ldclean::turner_effective_threads(threads);
        const auto start = std::chrono::steady_clock::now();
        double objective;
        if (profile) objective = ldclean::turner_min_cost_profiled(dfa, 100 * lambda, kernel, counters, 3, threads);
        else if (mode == "sparse") objective = ldclean::turner_min_cost(dfa, 100 * lambda, 3, threads);
        else if (mode == "dense_scalar") objective = ldclean::turner_min_cost_dense_right_normal(dfa, 100 * lambda, 3, threads);
        else objective = ldclean::turner_min_cost_dense_wavefront(dfa, 100 * lambda, 3, threads);
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        if (!std::isfinite(objective)) throw std::runtime_error("nonfinite objective");
        std::cout << std::setprecision(17) << "mode=" << mode << " lambda=" << lambda
                  << " threads_requested=" << threads << " threads_effective=" << effective
                  << " protein_aa=" << sequence.size() << " rna_nt=" << dfa.rna_length()
                  << " lattice_nodes=" << dfa.num_nodes() << " objective_units=" << objective
                  << " solve_seconds=" << seconds << " profiled=" << profile;
        if (profile) std::cout << " construction_seconds=" << counters.construction_seconds
                              << " closed_seconds=" << counters.closed_seconds
                              << " multiloop_seconds=" << counters.multiloop_seconds
                              << " external_seconds=" << counters.external_seconds
                              << " split_visits=" << counters.split_visits
                              << " split_eligible=" << counters.split_eligible;
        std::cout << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
