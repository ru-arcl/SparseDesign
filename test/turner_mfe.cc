// turner_mfe.cc — tiny driver: print the clean Turner-d0 design MFE (kcal/mol) for a protein.
// Usage: turner_mfe PROTEIN [lambda]   (run from lineardesign-clean/ for the data path)
#include <cstdio>
#include <cstdlib>

#include "../src/codon_table.h"
#include "../src/dfa.h"
#include "../src/fold_turner.h"

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: turner_mfe PROTEIN [lambda]\n"); return 1; }
    double lam = (argc > 2) ? atof(argv[2]) : 0.0;
    ldclean::CodonTable t("data/codon_usage_freq_table_human.csv");
    ldclean::DFA d(argv[1], t);
    double cost = ldclean::turner_min_cost(d, lam);
    printf("%.2f\n", cost / 100.0);
    return 0;
}
