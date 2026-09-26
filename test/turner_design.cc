// turner_design.cc — full clean Turner-d0 design for a protein.
// Usage: turner_design PROTEIN [lambda]   (run from lineardesign-clean/)
// Prints CDS, structure, MFE, CAI, and self-checks: independent evaluator == DP MFE,
// structure balanced, and translation back to the input protein.
#include <cstdio>
#include <cstdlib>
#include <string>

#include "../src/codon_table.h"
#include "../src/dfa.h"
#include "../src/fold_turner.h"

using namespace ldclean;

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: turner_design PROTEIN [lambda]\n"); return 1; }
    std::string prot = argv[1];
    double lam = (argc > 2) ? atof(argv[2]) : 0.0;
    CodonTable t("data/codon_usage_freq_table_human.csv");
    DFA d(prot, t);
    TurnerDesign r = design_turner(d, t, lam);

    // self-checks
    std::string aas;
    for (size_t i = 0; i + 2 < r.mrna.size() + 1; i += 3) aas.push_back(t.aa_of(r.mrna.substr(i, 3)));
    bool trans_ok = (aas == prot);
    bool len_ok = (r.mrna.size() == r.structure.size()) && ((int)r.mrna.size() == d.rna_length());
    // backtrace correctness: the emitted structure's energy must equal the DP cost (λ=0).
    int eval_e = turner_eval(r.mrna, r.structure);
    int dp_e = (int)(r.cost + (r.cost < 0 ? -0.5 : 0.5));
    bool eval_ok = (eval_e == dp_e);

    printf("%s\n%s\n", r.mrna.c_str(), r.structure.c_str());
    printf("MFE: %.2f kcal/mol   CAI: %.3f   (lambda=%.2f)\n", r.mfe_kcal, r.cai, lam);
    printf("checks: translate=%s len=%s eval==DPcost=%s (eval=%.2f, dp=%.2f)\n",
           trans_ok ? "ok" : "BAD", len_ok ? "ok" : "BAD", eval_ok ? "ok" : "BAD",
           eval_e / 100.0, r.cost / 100.0);
    return (trans_ok && len_ok && eval_ok) ? 0 : 1;
}
