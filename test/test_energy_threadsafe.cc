// Exercise the public energy API from many std::threads, including concurrent first access to
// the special-hairpin tables.  This is intentionally independent of OpenMP so it can also be
// run under ThreadSanitizer with just -pthread.
#include <atomic>
#include <cstdio>
#include <string>
#include <thread>
#include <vector>

#include "../src/energy.h"

namespace {

struct ExpectedHairpin {
    const char* span;
    int energy;
};

const ExpectedHairpin kTri[] = {{"CAACG", 680}, {"GUUAC", 690}};
const ExpectedHairpin kTet[] = {
    {"CAACGG", 550}, {"CCAAGG", 330}, {"CCACGG", 370}, {"CCCAGG", 340},
    {"CCGAGG", 350}, {"CCGCGG", 360}, {"CCUAGG", 370}, {"CCUCGG", 250},
    {"CUAAGG", 360}, {"CUACGG", 280}, {"CUCAGG", 370}, {"CUCCGG", 270},
    {"CUGCGG", 280}, {"CUUAGG", 350}, {"CUUCGG", 370}, {"CUUUGG", 370},
};
const ExpectedHairpin kHex[] = {
    {"ACAGUACU", 280}, {"ACAGUGAU", 360}, {"ACAGUGCU", 290}, {"ACAGUGUU", 180},
};

int nucleotide(char c) {
    if (c == 'A') return 1;
    if (c == 'C') return 2;
    if (c == 'G') return 3;
    if (c == 'U') return 4;
    return 0;
}

bool check_pairs() {
    static const int expected[5][5] = {
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 5},
        {0, 0, 0, 1, 0},
        {0, 0, 2, 0, 3},
        {0, 6, 0, 4, 0},
    };
    for (int a = -2; a <= 6; ++a) {
        for (int b = -2; b <= 6; ++b) {
            const int want = (a >= 0 && a <= 4 && b >= 0 && b <= 4) ? expected[a][b] : 0;
            if (ldclean::energy::pair_type(a, b) != want) return false;
            if (ldclean::energy::pairable(a, b) != (want != 0)) return false;
        }
    }
    return true;
}

template <size_t N>
bool check_hairpins(int size, const ExpectedHairpin (&expected)[N]) {
    const auto& got = ldclean::energy::special_hairpins(size);
    if (got.size() != N) return false;
    for (size_t i = 0; i < N; ++i) {
        const std::string span(expected[i].span);
        if (got[i].first != span || got[i].second != expected[i].energy) return false;
        if (ldclean::energy::special_index(span) != static_cast<int>(i)) return false;
        if (ldclean::energy::hairpin(size, nucleotide(span.front()), 1, 1,
                                     nucleotide(span.back()), static_cast<int>(i)) !=
            expected[i].energy)
            return false;
    }
    return true;
}

bool check_hairpins(int size) {
    if (size == 3) return check_hairpins(size, kTri);
    if (size == 4) return check_hairpins(size, kTet);
    if (size == 6) return check_hairpins(size, kHex);
    return ldclean::energy::special_hairpins(size).empty();
}

}  // namespace

int main() {
    constexpr int kThreadCount = 24;
    constexpr int kRounds = 250;
    const int sizes[] = {3, 4, 6, 5};
    std::atomic<int> ready(0);
    std::atomic<int> failures(0);
    std::atomic<bool> start(false);
    std::vector<std::thread> threads;
    threads.reserve(kThreadCount);

    for (int id = 0; id < kThreadCount; ++id) {
        threads.emplace_back([&, id] {
            ready.fetch_add(1, std::memory_order_release);
            while (!start.load(std::memory_order_acquire)) std::this_thread::yield();

            bool ok = true;
            for (int round = 0; round < kRounds && ok; ++round) {
                // init() remains in the API for compatibility and must be harmless even when
                // old callers invoke it repeatedly and concurrently.
                ldclean::energy::init();
                ok = check_pairs() && check_hairpins(sizes[(id + round) % 4]);
            }
            if (!ok) failures.fetch_add(1, std::memory_order_relaxed);
        });
    }

    while (ready.load(std::memory_order_acquire) != kThreadCount) std::this_thread::yield();
    start.store(true, std::memory_order_release);
    for (auto& thread : threads) thread.join();

    const int failed = failures.load(std::memory_order_relaxed);
    std::printf("energy thread-safety: %s (%d threads, %d rounds)\n",
                failed == 0 ? "PASS" : "FAIL", kThreadCount, kRounds);
    return failed == 0 ? 0 : 1;
}
