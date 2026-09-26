// test_codon.cc — sanity checks for codon_table against the original tool's reported values.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

#include "../src/codon_table.h"

static int failures = 0;
static void check_near(const char* name, double got, double want, double tol) {
    bool ok = std::fabs(got - want) <= tol;
    printf("%-28s got=%.4f want=%.4f  %s\n", name, got, want, ok ? "OK" : "FAIL");
    if (!ok) ++failures;
}

static bool invalid_frequency_rejected(const std::string& label, const std::string& frequency) {
    const std::string path = "/tmp/ldclean-test-codon-" + label + ".csv";
    {
        std::ofstream out(path.c_str());
        out << "AAA,K," << frequency << "\n";
    }
    bool rejected = false;
    try {
        ldclean::CodonTable invalid(path);
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    std::remove(path.c_str());
    printf("invalid frequency %-8s %s\n", label.c_str(), rejected ? "rejected" : "ACCEPTED");
    if (!rejected) ++failures;
    return rejected;
}

int main() {
    ldclean::CodonTable t("data/codon_usage_freq_table_human.csv");

    // The original tool reported CAI 0.734 for this P15421 design at lambda=0
    // (lineardesign/reference/baseline_so_outputs.txt). Includes the stop codon.
    std::string p15421_l0 =
        "AUGUAUGGCAAGAUCAUCUUCGUCCUGCUGCUCUCCGGGAUCGUGUCGAUCUCGGCGAGCAGCACGACGGGGGUGGCC"
        "AUGCAUACGAGUACCAGCAGUAGCGUGACUAAGAGUUAUAUAUCCUCACAGACCAACGGCAUCACCUUGAUAAAUUGGUGG"
        "GCGAUGGCCCGCGUAAUUUUCGAGGUGAUGCUGGUGGUCGUGGGGAUGAUAAUUCUUAUCAGCUACUGCAUUCGU";
    check_near("CAI(P15421, lambda=0)", t.cai(p15421_l0), 0.734, 0.0015);

    // A pure max-codon sequence should have CAI very close to 1.0.
    // MFLMVF at lambda=0 from the original was AUGUUCCUGAUGGUGUUC, CAI 1.000.
    check_near("CAI(MFLMVF design)", t.cai("AUGUUCCUGAUGGUGUUC"), 1.000, 0.0015);

    // w(c) in (0,1]; the most-frequent synonymous codon has w == 1.
    printf("\nspot-check w(c):\n");
    for (const auto& c : {std::string("AUG"), std::string("UUC"), std::string("CUG"),
                          std::string("GUG"), std::string("UAA")}) {
        printf("  w(%s) = %.3f  (aa=%c)\n", c.c_str(), t.w_of(c), t.aa_of(c));
    }
    // Leucine should have 6 synonymous codons; check counts.
    printf("\nsynonymous-codon counts: L=%zu R=%zu S=%zu M=%zu *=%zu\n",
           t.options_for('L').size(), t.options_for('R').size(), t.options_for('S').size(),
           t.options_for('M').size(), t.options_for('*').size());

    bool finite_weights = true;
    for (char aa : std::string("*ACDEFGHIKLMNPQRSTVWY"))
        for (const auto& option : t.options_for(aa))
            finite_weights = finite_weights && std::isfinite(option.freq) && option.freq >= 0.0 &&
                             std::isfinite(option.w) && option.w > 0.0 && option.w <= 1.0 &&
                             std::isfinite(option.neg_log_w) && option.neg_log_w >= 0.0;
    printf("all parsed costs finite/nonnegative: %s\n", finite_weights ? "OK" : "FAIL");
    if (!finite_weights) ++failures;

    invalid_frequency_rejected("nan", "nan");
    invalid_frequency_rejected("inf", "inf");
    invalid_frequency_rejected("negative", "-0.1");

    // A zero-frequency codon remains accepted and receives the documented finite floor.
    const std::string zero_path = "/tmp/ldclean-test-codon-zero.csv";
    {
        std::ofstream out(zero_path.c_str());
        out << "AAA,K,0\n";
    }
    bool zero_ok = false;
    try {
        ldclean::CodonTable zero(zero_path);
        const double w = zero.w_of("AAA");
        zero_ok = std::isfinite(w) && w > 0.0 && w <= 1.0 &&
                  std::isfinite(zero.options_for('K')[0].neg_log_w);
    } catch (const std::exception&) {
        zero_ok = false;
    }
    std::remove(zero_path.c_str());
    printf("zero frequency floor: %s\n", zero_ok ? "OK" : "FAIL");
    if (!zero_ok) ++failures;

    printf("\n%s\n", failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED");
    return failures == 0 ? 0 : 1;
}
