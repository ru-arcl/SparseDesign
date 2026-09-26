// Focused ablation check: the validation-only dense bottom-up wavefront must reproduce both the
// retained dense top-down oracle and the production sparse wavefront at every thread count.
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "../src/codon_table.h"
#include "../src/dfa.h"
#include "../src/fold_turner.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern int ML_BASE37;

using namespace ldclean;

struct TestCase {
    const char* protein;
    double lambda;
    int ml_base;
};

static int observe_team(int requested) {
    int observed = 1;
#ifdef _OPENMP
#pragma omp parallel num_threads(requested)
    {
#pragma omp single
        observed = omp_get_num_threads();
    }
#else
    (void)requested;
#endif
    return observed;
}

int main() {
    CodonTable table("data/codon_usage_freq_table_human.csv");
    const int saved_ml_base = ML_BASE37;
    const std::vector<TestCase> cases = {
        {"MFL", 0.0, 0},
        {"SSSS", 100.0, 0},
        {"MSWFWFM", 300.0, 0},
        // These production-level regression cases select a multiloop; the second also exercises
        // an all-unpaired multiloop suffix under a nonzero per-nucleotide penalty.
        {"WIPKKSAWCFEQQLCNLYSQLFMGV", 400.0, 0},
        {"WIPKKSAWCFEQQLCNLYSQGFMGV", 400.0, 37},
    };
    const std::vector<int> threads = {1, 2, 4};

    int failures = 0;
    // The solver policy cap is not necessarily the runtime team size: OMP_THREAD_LIMIT,
    // OMP_DYNAMIC, nested execution, and non-OpenMP builds may reduce it. Verify the public
    // reporting helper against an independent team using its reported request. In particular,
    // running this test with OMP_THREAD_LIMIT=1 guards against claiming j4 while running serially.
#ifdef _OPENMP
    const int saved_dynamic = omp_get_dynamic();
    omp_set_dynamic(0);  // make the two independent probes deterministic
#endif
    for (int requested : threads) {
        const int reported = turner_effective_threads(requested);
        const int observed = observe_team(reported);
        const bool ok = reported == observed;
        std::printf("thread-report requested=%d reported=%d observed=%d %s\n",
                    requested, reported, observed, ok ? "OK" : "FAIL");
        if (!ok) ++failures;
    }
#ifdef _OPENMP
    omp_set_dynamic(saved_dynamic);
#endif
    for (const TestCase& test : cases) {
        ML_BASE37 = test.ml_base;
        DFA lattice(test.protein, table);
        const double reference = turner_min_cost_dense_reference(lattice, test.lambda, 3);
        const double sparse = turner_min_cost(lattice, test.lambda, 3, 1);
        bool ok = std::fabs(reference - sparse) < 1e-6;
        std::printf("protein=%s lambda=%.2f ML_BASE37=%d dense_ref=%.9f sparse=%.9f",
                    test.protein, test.lambda, test.ml_base, reference, sparse);
        for (int nt : threads) {
            const double wave = turner_min_cost_dense_wavefront(lattice, test.lambda, 3, nt);
            const double right = turner_min_cost_dense_right_normal(lattice, test.lambda, 3, nt);
            if (std::fabs(reference - right) >= 1e-6) ok = false;
            std::printf(" dense_wave_j%d=%.9f", nt, wave);
            if (std::fabs(reference - wave) >= 1e-6) ok = false;
        }
        std::printf(" %s\n", ok ? "OK" : "FAIL");
        if (!ok) ++failures;
    }
    // Product-lattice regression: several automaton states share a nucleotide layer but have no
    // epsilon path between them, exercising the dense kernel's same-layer endpoint handling.
    ML_BASE37 = 0;
    DFA constrained("SSSS", table, {"UCUUCU", "AGCAGC"});
    const double constrained_reference = turner_min_cost_dense_reference(constrained, 100.0, 3);
    const double constrained_sparse = turner_min_cost(constrained, 100.0, 3, 1);
    bool constrained_ok = std::fabs(constrained_reference - constrained_sparse) < 1e-6;
    std::printf("protein=SSSS constrained=yes lambda=100.00 dense_ref=%.9f sparse=%.9f",
                constrained_reference, constrained_sparse);
    for (int nt : threads) {
        const double wave = turner_min_cost_dense_wavefront(constrained, 100.0, 3, nt);
        const double right = turner_min_cost_dense_right_normal(constrained, 100.0, 3, nt);
        if (std::fabs(constrained_reference - right) >= 1e-6) constrained_ok = false;
        std::printf(" dense_wave_j%d=%.9f", nt, wave);
        if (std::fabs(constrained_reference - wave) >= 1e-6) constrained_ok = false;
    }
    std::printf(" %s\n", constrained_ok ? "OK" : "FAIL");
    if (!constrained_ok) ++failures;

    ML_BASE37 = saved_ml_base;
    std::printf("dense-wavefront equivalence: %zu cases, %d failures\n",
                cases.size() + 1, failures);
    return failures == 0 ? 0 : 1;
}
