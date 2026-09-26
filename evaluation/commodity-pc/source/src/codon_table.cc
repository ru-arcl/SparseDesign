// codon_table.cc — see codon_table.h.
#include "codon_table.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace ldclean {

namespace {
// Normalize a raw codon field to uppercase RNA (T -> U).
std::string normalize_codon(std::string c) {
    for (char& ch : c) {
        ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        if (ch == 'T') ch = 'U';
    }
    return c;
}
std::string strip(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}
}  // namespace

CodonTable::CodonTable(const std::string& csv_path) {
    std::ifstream in(csv_path);
    if (!in) throw std::runtime_error("cannot open codon usage table: " + csv_path);

    // First pass: read (codon, aa, freq); track max synonymous freq per aa.
    std::map<char, double> max_freq;
    std::vector<std::tuple<std::string, char, double>> rows;
    std::string line;
    while (std::getline(in, line)) {
        // strip UTF-8 BOM if present on the first line
        if (line.size() >= 3 && (unsigned char)line[0] == 0xEF &&
            (unsigned char)line[1] == 0xBB && (unsigned char)line[2] == 0xBF)
            line = line.substr(3);
        line = strip(line);
        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);
        std::string codon_f, aa_f, freq_f;
        if (!std::getline(ss, codon_f, ',')) continue;
        if (!std::getline(ss, aa_f, ',')) continue;
        if (!std::getline(ss, freq_f, ',')) continue;

        std::string codon = normalize_codon(strip(codon_f));
        std::string aa_s = strip(aa_f);
        if (codon.size() != 3 || aa_s.empty()) continue;
        char aa = aa_s[0];
        const std::string frequency = strip(freq_f);
        std::size_t consumed = 0;
        double freq = 0.0;
        try {
            freq = std::stod(frequency, &consumed);
        } catch (const std::exception&) {
            throw std::runtime_error("invalid codon frequency for " + codon + ": '" +
                                     frequency + "'");
        }
        // Zero remains a supported representation for an unused codon and is floored to 1e-9
        // below. NaN/Inf/negative values would otherwise leak into DFA edge weights and make the
        // folding DP's ordering and sentinel invariants undefined.
        if (consumed != frequency.size() || !std::isfinite(freq) || freq < 0.0)
            throw std::runtime_error("invalid codon frequency for " + codon + ": '" +
                                     frequency + "'");

        rows.emplace_back(codon, aa, freq);
        auto it = max_freq.find(aa);
        if (it == max_freq.end() || freq > it->second) max_freq[aa] = freq;
    }
    if (rows.empty()) throw std::runtime_error("no codon rows parsed from: " + csv_path);

    // Second pass: compute w(c) = freq / max synonymous freq.
    for (auto& row : rows) {
        const std::string& codon = std::get<0>(row);
        char aa = std::get<1>(row);
        double freq = std::get<2>(row);
        double fmax = max_freq[aa];
        double w = (fmax > 0.0) ? (freq / fmax) : 0.0;
        if (w <= 0.0) w = 1e-9;  // guard against log(0) for never-used codons

        CodonOption opt{codon, freq, w, -std::log(w)};
        aa_to_options_[aa].push_back(opt);
        codon_to_aa_[codon] = aa;
        codon_to_w_[codon] = w;
    }
    // Stable ordering of options (by codon) for determinism.
    for (auto& kv : aa_to_options_)
        std::sort(kv.second.begin(), kv.second.end(),
                  [](const CodonOption& a, const CodonOption& b) { return a.codon < b.codon; });
}

const std::vector<CodonOption>& CodonTable::options_for(char aa) const {
    auto it = aa_to_options_.find(aa);
    if (it == aa_to_options_.end())
        throw std::runtime_error(std::string("unknown amino acid: ") + aa);
    return it->second;
}

char CodonTable::aa_of(const std::string& codon) const {
    auto it = codon_to_aa_.find(codon);
    return it == codon_to_aa_.end() ? '\0' : it->second;
}

double CodonTable::w_of(const std::string& codon) const {
    auto it = codon_to_w_.find(codon);
    return it == codon_to_w_.end() ? 0.0 : it->second;
}

double CodonTable::cai(const std::string& rna) const {
    if (rna.empty() || rna.size() % 3 != 0) return 0.0;
    int ncodons = static_cast<int>(rna.size() / 3);
    double sum_log = 0.0;
    for (int i = 0; i < ncodons; ++i) {
        std::string c = rna.substr(i * 3, 3);
        double w = w_of(c);
        if (w <= 0.0) w = 1e-9;
        sum_log += std::log(w);
    }
    return std::exp(sum_log / ncodons);
}

}  // namespace ldclean
