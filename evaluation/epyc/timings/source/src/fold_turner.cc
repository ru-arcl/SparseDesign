// fold_turner.cc — see fold_turner.h. Turner-d0 lattice DP with backtrace + evaluator.
#include "fold_turner.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <sys/stat.h>
#include <thread>

#include "energy.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// Layout selection. The production Makefile defines LDCLEAN_SQUARE_MEMOS because the full n*n
// arena is fastest at both one thread and 16 threads. Omitting it selects the forward-only packed
// arena used by the explicitly named low-memory build. LDCLEAN_NO_SAME_LAYER_PRESEED and
// LDCLEAN_PACKED_RANK_MAP retain historical validation ablations. The lazy same-layer mode is
// intentionally unavailable for a packed arena, where a reverse-layer probe has no storage, or
// with forward-key assertions enabled.
#if defined(LDCLEAN_NO_SAME_LAYER_PRESEED) && !defined(LDCLEAN_SQUARE_MEMOS)
#error "LDCLEAN_NO_SAME_LAYER_PRESEED requires LDCLEAN_SQUARE_MEMOS"
#endif
#if defined(LDCLEAN_NO_SAME_LAYER_PRESEED) && defined(LDCLEAN_VALIDATE_FORWARD_KEYS)
#error "LDCLEAN_NO_SAME_LAYER_PRESEED is incompatible with LDCLEAN_VALIDATE_FORWARD_KEYS"
#endif
#if defined(LDCLEAN_PACKED_RANK_MAP) && defined(LDCLEAN_SQUARE_MEMOS)
#error "LDCLEAN_PACKED_RANK_MAP is incompatible with LDCLEAN_SQUARE_MEMOS"
#endif

namespace ldclean {
namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();
// The most-negative finite double is below every possible Turner objective: lambda and all DFA
// weights are nonnegative, while the RNA length is int-bounded and every energy term is an int.
// Keeping an ordinary comparison here avoids an extra computed-flag load in the hottest memos.
constexpr double kUncomputed = std::numeric_limits<double>::lowest();
constexpr int kMaxLoop = 30;  // Vienna internal/bulge loop cap
// Ceiling on effective parallelism, detected from the host at runtime (was a hardcoded 16 tuned
// for one 2-socket box). This DP scales near-linearly within a single CPU socket, but threads
// spanning two NUMA sockets regress it (remote-memory traffic on the shared memo tables). So:
// use ALL cores on a single-socket host; cap at one node's worth of cores on a multi-socket host.
inline int useful_thread_ceiling() {
    unsigned hc = std::thread::hardware_concurrency();
    if (hc == 0) hc = 1;
    // Count NUMA nodes via /sys/devices/system/node/nodeN (Linux); >1 => multi-socket.
    int nodes = 0;
    for (int i = 0; i < 256; ++i) {
        std::string p = "/sys/devices/system/node/node" + std::to_string(i);
        struct stat st;
        if (::stat(p.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) ++nodes; else break;
    }
    if (nodes > 1) return std::max(1, (int)hc / nodes);
    return (int)hc;
}

struct PairCost { int nt5; int nt3; double cost; };
struct GapOpt { int first; int last; double cai; int x2; int y2; };  // unpaired stretch option
// One dense multiloop branch used to replay the original recurrence during backtrace: start-node a
// pairs to split-node c with branch cost `bcost` = closed-cost + multi_branch(nt5,nt3).
struct BranchEntry { int c; int nt5; int nt3; double bcost; };
// Exact sparse-multiloop candidate [start,end]. For a fixed end node, only a direct branch whose
// cost is strictly better than every partitionable/unpaired alternative can be needed as the
// rightmost branch of m2. At most one (the cheapest pair type) is retained per node pair.
struct MultiCandidate { int start; double bcost; };

// Compact descriptions used while reconstructing the winning path after the cost DP. Unpaired
// stretches are stored as the chosen GapOpt + its span (x<0 = empty side), then reconstructed to
// nucleotides only once the winning decision is known.
struct TLDetail {
    int ia, ib, ntp, ntq;          // inner pair (nodes + closing nts)
    int x5, y5, x3, y3;            // 5'/3' unpaired spans (x<0 = no unpaired nts that side)
    GapOpt g5, g3;                 // chosen gap option per side (reconstructed via gap_string)
};
struct HpDec {                     // hairpin loop, deferred
    bool special;                  // special tri/tetra/hexaloop -> use `loop`; else gap_string(ap,bp,g)
    std::string loop;              // special loop (<= 6 nt; built cheaply by hp_enum, no mincai_path)
    GapOpt g;                      // generic loop option
};
struct ClosedDec { int ap, bp, kind; HpDec hp; TLDetail tl; };   // kind 0=hp 1=two 2=multi
struct NodeDec { int kind; int nt, to; int c, nt5, nt3; };       // ext/m1/m2 choice

inline bool eq(double a, double b) { return std::fabs(a - b) < 1e-4; }

class TurnerSolver {
public:
    TurnerSolver(const DFA& d, double lambda, int min_loop, int threads)
        : d_(d), lambda_(lambda), min_loop_(min_loop),
          num_threads_(std::max(1, std::min(threads, useful_thread_ceiling()))) {
        if (!std::isfinite(lambda_) || lambda_ < 0.0)
            throw std::invalid_argument("Turner lambda must be finite and nonnegative");
        if (min_loop_ < 0)
            throw std::invalid_argument("Turner minimum loop size must be nonnegative");
        energy::init();
        // Every memoized interval is layer-forward: pos(a) <= pos(b). Stable DFA construction
        // makes node IDs position-major without changing within-layer or edge iteration order.
        // The packed low-memory arena stores each start-node row from the beginning of a's layer
        // onward, including the full same-layer Cartesian block needed by the base cases below.
        // The speed-oriented production build defines LDCLEAN_SQUARE_MEMOS instead.
        n_ = d_.num_nodes();
#ifdef LDCLEAN_SQUARE_MEMOS
        const size_t cells = static_cast<size_t>(n_) * static_cast<size_t>(n_);
#else
        row_base_.resize(n_);
#ifdef LDCLEAN_PACKED_RANK_MAP
        rank_by_position_.assign(n_, -1);
#endif
        size_t next_node = 0;
        size_t layer_square_cells = 0;
        for (int p = 0; p <= d_.rna_length(); ++p) {
            const size_t width = d_.nodes_at(p).size();
            layer_square_cells += width * width;
            for (int node : d_.nodes_at(p)) {
                if (node < 0 || node >= n_ || node != static_cast<int>(next_node) ||
                    d_.pos_of(node) != p)
                    throw std::logic_error("DFA node IDs are not position-major");
#ifdef LDCLEAN_PACKED_RANK_MAP
                rank_by_position_[node] = static_cast<int>(next_node);
#endif
                ++next_node;
            }
        }
        if (next_node != static_cast<size_t>(n_))
            throw std::logic_error("DFA layers do not cover every node ID");

        size_t cells = 0;
        size_t layer_first = 0;
        for (int p = 0; p <= d_.rna_length(); ++p) {
            const size_t width = d_.nodes_at(p).size();
            const size_t row_length = static_cast<size_t>(n_) - layer_first;
            for (int node : d_.nodes_at(p)) {
                // key(node,b) = row_offset[node] + b - layer_first. Folding the layer_first
                // subtraction into this adjusted base leaves one load and one add in the hot
                // lookup.
                row_base_[node] = cells - layer_first;
                cells += row_length;
            }
            layer_first += width;
        }
        const size_t expected_cells =
            (static_cast<size_t>(n_) * n_ + layer_square_cells) / 2;
        if (layer_first != static_cast<size_t>(n_) || cells != expected_cells)
            throw std::logic_error("invalid forward-pair memo layout");
#endif
        memo_cells_ = cells;

        // A sentinel in each scalar value replaces both per-lookup hashing and a separate
        // computed flag. ext remains 1-D because its right endpoint is invariant.
        mincai_v_.assign(cells, kUncomputed);
        // Every external-loop recursion keeps b fixed at d_.final(), so ext is genuinely 1-D.
        ext_v_.assign(n_, kUncomputed);
        m1_v_.assign(cells, kUncomputed);
        m2_v_.assign(cells, kUncomputed);
        closed_v_.resize(cells); closed_c_.assign(cells, 0);
        gap_v_.resize(cells);    gap_c_.assign(cells, 0);
        // Every edge advances by exactly one nucleotide layer. Therefore distinct states on the
        // same layer have no path between them: their all-unpaired fragment is infeasible, as are
        // all same-layer one/multiple-branch multiloop fragments. Seed those base cases once for
        // every solver mode. Besides making the layer-order invariant explicit, this prevents an
        // unreachable top-down probe from walking past its fixed right endpoint.
#ifndef LDCLEAN_NO_SAME_LAYER_PRESEED
        for (int p = 0; p <= d_.rna_length(); ++p)
            for (int a : d_.nodes_at(p))
                for (int b : d_.nodes_at(p)) {
                    const size_t k = key(a, b);
                    mincai_v_[k] = (a == b) ? 0.0 : kInf;
                    m1_v_[k] = kInf;
                    m2_v_[k] = kInf;
                }
#endif
        // Backtrace replays the original dense recurrence to preserve its exact tie order. Build
        // those few required branch rows on demand; the forward DP uses mcand_end_ instead.
        brow_.assign(n_, {});
        brow_hi_.assign(n_, -1);
        mcand_end_.assign(n_, {});
    }

    double solve() {
        // The sparse wavefront removes the cubic branch rescan that previously made a one-thread
        // bottom-up fill slower than lazy recursion. Use it for both serial and parallel solves;
        // each diagonal remains independently parallel when num_threads_ > 1.
        fill_bottom_up<true, false>();
        return ext(d_.start(), d_.final());
    }

    // Validation-only oracle: without the bottom-up sparse fill, closed() reaches the retained
    // top-down m1()/m2() functions, which enumerate every direct multiloop branch.
    double solve_dense_reference() { return ext(d_.start(), d_.final()); }

    // Validation/ablation path: the dense recurrence on the same compact wavefront infrastructure
    // as solve(), differing only in its eager all-branch rows and fused dense multiloop kernel.
    double solve_dense_wavefront() {
        fill_bottom_up_dense<false>();
        return ext(d_.start(), d_.final());
    }

    double solve_dense_right_normal() {
        fill_bottom_up<false, false>();
        return ext(d_.start(), d_.final());
    }

    double solve_profiled(TurnerKernel kernel, TurnerProfile& out) {
        split_visits_.assign(n_, 0);
        split_eligible_.assign(n_, 0);
        if (kernel == TurnerKernel::Sparse) fill_bottom_up<true, true>();
        else if (kernel == TurnerKernel::DenseRightNormal) fill_bottom_up<false, true>();
        else fill_bottom_up_dense<true>();
        out.closed_seconds = profile_.closed_seconds;
        out.multiloop_seconds = profile_.multiloop_seconds;
        out.external_seconds = profile_.external_seconds;
        for (int i = 0; i < n_; ++i) {
            out.split_visits += split_visits_[i];
            out.split_eligible += split_eligible_[i];
        }
        return ext(d_.start(), d_.final());
    }

#ifdef LDCLEAN_TEST_CELL_PARITY
    // Test-only exposure of completed memo cells; absent from production/API builds.
    double test_m1(int a, int b) { return m1(a, b); }
    double test_m2(int a, int b) { return m2(a, b); }
    const std::vector<PairCost>& test_closed(int a, int b) { return closed(a, b); }
#endif

    void sparse_stats(TurnerSparseStats& out) const {
        out = TurnerSparseStats{};
        out.rna_length = d_.rna_length();
        out.lattice_nodes = static_cast<size_t>(n_);
        const int L = d_.rna_length();
        for (int p = 0; p < L; ++p)
            for (int q = p + 1; q <= L; ++q)
                out.interval_cells += d_.nodes_at(p).size() * d_.nodes_at(q).size();
        // Reverse cells in the square arena remain empty; the packed arena omits them. A flat scan
        // therefore counts feasible direct intervals in either layout without another DP pass.
        for (size_t k = 0; k < closed_v_.size(); ++k)
            if (!closed_v_[k].empty()) ++out.direct_intervals;
        for (const auto& by_end : mcand_end_) {
            out.candidates += by_end.size();
            out.candidate_capacity += by_end.capacity();
            out.max_candidates_per_end = std::max(out.max_candidates_per_end, by_end.size());
        }
        out.candidate_storage_bytes = mcand_end_.size() * sizeof(mcand_end_[0]) +
                                      out.candidate_capacity * sizeof(MultiCandidate);
    }

    void backtrace(std::string& seq, std::string& st) { bt_ext(d_.start(), d_.final(), seq, st); }

private:
    // Flat memo index. The default position-major packed path needs one offset load and one add;
    // LDCLEAN_SQUARE_MEMOS restores the historical multiply/add for validation A/Bs.
    size_t key(int a, int b) const {
#ifdef LDCLEAN_VALIDATE_FORWARD_KEYS
        assert(a >= 0 && a < n_ && b >= 0 && b < n_);
        assert(d_.pos_of(a) <= d_.pos_of(b));
#ifdef LDCLEAN_PACKED_RANK_MAP
        assert(rank_by_position_[b] >= 0);
#endif
#endif
#ifdef LDCLEAN_SQUARE_MEMOS
        const size_t k = static_cast<size_t>(a) * static_cast<size_t>(n_) +
                         static_cast<size_t>(b);
#elif defined(LDCLEAN_PACKED_RANK_MAP)
        const size_t k = row_base_[a] + static_cast<size_t>(rank_by_position_[b]);
#else
        const size_t k = row_base_[a] + static_cast<size_t>(b);
#endif
#ifdef LDCLEAN_VALIDATE_FORWARD_KEYS
        assert(k < memo_cells_);
#endif
        return k;
    }

    // Bottom-up evaluation: every cell (a,b) of span-gap g depends only on cells of gap < g, so
    // filling in increasing-gap order lets each memoized function compute exactly its own cell
    // (all sub-cells already cached) with no deep recursion -- better cache locality than lazy
    // top-down, and embarrassingly parallel within each gap.
    template <bool Prune, bool Profile>
    void fill_bottom_up() {
        const int L = d_.rna_length();
        const int nt = num_threads_;
        // The bottom-up path uses exact sparse-multiloop candidate lists keyed by right endpoint.
        // Candidates are produced in phase 2 after every dependency at a smaller gap is complete.
        // For a fixed gap, a right endpoint b belongs to only one p iteration, so appends to
        // mcand_end_[b] are lock-free. brow_ is reserved for on-demand backtrace only.
        for (int g = 1; g <= L; ++g) {
            auto phase_start = Profile ? std::chrono::steady_clock::now() :
                                         std::chrono::steady_clock::time_point{};
            // The three phases are barriered (closed -> m1/m2 -> ext) and all cells of a phase are
            // independent (they read only already-filled smaller-or-equal gaps), so each phase
            // parallelizes over the start position p; memoized sub-calls hit the filled tables.
            const int pmax = L - g;
            // Phase 1: mincai/gap_options/closed -- depend only on gaps < g.
#pragma omp parallel for schedule(dynamic) num_threads(nt)
            for (int p = 0; p <= pmax; ++p)
                for (int a : d_.nodes_at(p))
                    for (int b : d_.nodes_at(p + g)) {
                        mincai(a, b); gap_options(a, b); closed(a, b);
                    }
            if (Profile) record_phase(profile_.closed_seconds, phase_start);
            // Phase 2: m1/m2 for every gap-g cell (needed generally, as inner spans of enclosing
            // multiloops). ext is handled separately below.
#pragma omp parallel for schedule(dynamic) num_threads(nt)
            for (int p = 0; p <= pmax; ++p)
                for (int a : d_.nodes_at(p))
                    for (int b : d_.nodes_at(p + g)) m1m2_right_normal<Prune, Profile>(a, b);
            if (Profile) record_phase(profile_.multiloop_seconds, phase_start);
            // Phase 3: ext -- every ext recursion keeps its second arg fixed at final(), and so
            // does backtrace, so ext(a,b) is *read* only when b==final(). Fill just that column:
            // the gap-g cells with b==final are exactly the nodes at position (L-g). This turns
            // the ext phase from O(n^2) filled cells into O(n), removing most of the bottom-up
            // fill's excess work over the (reachable-only) lazy top-down path.
#pragma omp parallel for schedule(dynamic) num_threads(nt)
            for (int ai = 0; ai < (int)d_.nodes_at(L - g).size(); ++ai)
                ext(d_.nodes_at(L - g)[ai], d_.final());
            if (Profile) record_phase(profile_.external_seconds, phase_start);
        }
    }

    // Dense-wavefront ablation used only by turner_min_cost_dense_wavefront(). Phase 1 eagerly
    // builds each start node's complete direct-branch row in the exact position/node/pair order of
    // each_branch(); phase 2 scans that row for every interval instead of pruning candidates.
    // All scalar/vector memo layouts and the external-loop phase are shared with the sparse path.
    template <bool Profile>
    void fill_bottom_up_dense() {
        const int L = d_.rna_length();
        const int nt = num_threads_;
#ifdef LDCLEAN_NO_SAME_LAYER_PRESEED
        // Historical square-memo ablation A left same-layer cells lazy. The dense fused kernel
        // nevertheless requires its empty-span m1/m2 dependencies to be explicit before phase 1.
        for (int p = 0; p <= L; ++p)
            for (int a : d_.nodes_at(p))
                for (int b : d_.nodes_at(p)) {
                    m1_v_[key(a, b)] = kInf;
                    m2_v_[key(a, b)] = kInf;
                }
#endif
        for (int g = 1; g <= L; ++g) {
            auto phase_start = Profile ? std::chrono::steady_clock::now() :
                                         std::chrono::steady_clock::time_point{};
            const int pmax = L - g;
#pragma omp parallel for schedule(dynamic) num_threads(nt)
            for (int p = 0; p <= pmax; ++p)
                for (int a : d_.nodes_at(p))
                    for (int b : d_.nodes_at(p + g)) {
                        mincai(a, b);
                        gap_options(a, b);
                        for (const PairCost& pc : closed(a, b))
                            brow_[a].push_back(BranchEntry{
                                b, pc.nt5, pc.nt3,
                                pc.cost + energy::multi_branch(pc.nt5, pc.nt3)});
                    }
            if (Profile) record_phase(profile_.closed_seconds, phase_start);
#pragma omp parallel for schedule(dynamic) num_threads(nt)
            for (int p = 0; p <= pmax; ++p)
                for (int a : d_.nodes_at(p))
                    for (int b : d_.nodes_at(p + g)) m1m2_dense<Profile>(a, b);
            if (Profile) record_phase(profile_.multiloop_seconds, phase_start);
#pragma omp parallel for schedule(dynamic) num_threads(nt)
            for (int ai = 0; ai < static_cast<int>(d_.nodes_at(L - g).size()); ++ai)
                ext(d_.nodes_at(L - g)[ai], d_.final());
            if (Profile) record_phase(profile_.external_seconds, phase_start);
        }
        // The rows are now complete. If decision replay is added to this validation path later,
        // prevent each_branch() from treating an eager row as an incomplete lazy row.
        std::fill(brow_hi_.begin(), brow_hi_.end(), L);
    }

    // ---------------- cost DP ----------------
    double mincai(int a, int b) {
        if (a == b) return 0.0;
        size_t k = key(a, b);
        if (mincai_v_[k] != kUncomputed) return mincai_v_[k];
        double best = kInf;
        for (const Edge& e : d_.out_edges(a)) {
            double sub = mincai(e.to, b);
            if (sub < kInf) best = std::min(best, lambda_ * e.weight + sub);
        }
        mincai_v_[k] = best;
        return best;
    }

    // Minimum codon cost of an all-unpaired multiloop stretch, including the Turner per-base
    // term. Keeping this separate from mincai() prevents the dense recurrence and its replayed
    // backtrace decision from silently diverging when that term is nonzero.
    double multi_unpaired_gap(int a, int b) {
        double cost = mincai(a, b);
        if (cost < kInf)
            cost += static_cast<double>(d_.pos_of(b) - d_.pos_of(a)) *
                    energy::multi_unpaired();
        return cost;
    }

    const std::vector<GapOpt>& gap_options(int x, int y) {
        size_t k = key(x, y);
        if (gap_c_[k]) return gap_v_[k];
        std::vector<GapOpt>& opts = gap_v_[k];
        int dp = d_.pos_of(y) - d_.pos_of(x);
        if (dp == 0) {
            // A zero-length DFA path exists only from a state to itself. Product-lattice states
            // on the same nucleotide layer are distinct histories, not epsilon-connected aliases.
            if (x == y) opts.push_back({-1, -1, 0.0, x, y});
        } else if (dp == 1) {
            for (const Edge& e : d_.out_edges(x))
                if (e.to == y) opts.push_back({e.nuc, e.nuc, lambda_ * e.weight, y, x});
        } else {
            for (const Edge& ef : d_.out_edges(x))
                for (const Edge& el : d_.in_edges(y)) {
                    double mid = mincai(ef.to, el.from);
                    if (mid >= kInf) continue;
                    opts.push_back({ef.nuc, el.nuc, lambda_ * (ef.weight + el.weight) + mid,
                                    ef.to, el.from});
                }
        }
        gap_c_[k] = 1;
        return opts;
    }


    // hairpin closed by (nt5@i, nt3@j-1); inner span ap..bp. If best_loop, store argmin loop nts.
    double hairpin_cost(int ap, int bp, int nt5, int nt3, HpDec* dec) {
        int size = d_.pos_of(bp) - d_.pos_of(ap);
        if (size < min_loop_) return kInf;
        // Sizes 3/4/6 have special hairpins (tri/tetra/hexaloops) whose energy depends on the
        // FULL loop sequence. The generic "loop-independent energy + min-CAI middle" decomposition
        // is therefore UNSOUND for them: it can select a min-CAI loop that is in fact a special
        // hairpin and score it with the (cheaper) generic energy, while the structure evaluator
        // charges the special energy -> a silently suboptimal design. These sizes are tiny
        // (<= 4^6 DFA paths), so enumerate them explicitly with the true (special-aware) energy,
        // exactly matching turner_eval(). (Larger generic loops use the fast decomposition below,
        // where no special hairpins exist.)
        if (size == 3 || size == 4 || size == 6) {
            double best = kInf;
            std::string loop, bestloop;
            hp_enum(ap, bp, size, nt5, nt3, loop, 0.0, best, dec ? &bestloop : nullptr);
            if (dec) { dec->special = true; dec->loop = std::move(bestloop); }
            return best;
        }
        // Other sizes: no special hairpins. Generic energy depends only on the first/last loop nt
        // (terminal mismatch), so minimize CAI over the middle via gap_options. Defer the loop
        // string: store the winning GapOpt and rebuild it at backtrace (no mincai_path here).
        double best = kInf;
        for (auto& g : gap_options(ap, bp)) {
            int e = energy::hairpin(size, nt5, g.first, g.last, nt3, -1);
            double v = (double)e + g.cai;
            if (v < best) { best = v; if (dec) { dec->special = false; dec->g = g; } }
        }
        return best;
    }

    // Enumerate every loop sequence of a small special-hairpin size (3/4/6) reachable through the
    // DFA from `cur` toward `bp`, scoring each with its TRUE (special-aware) hairpin energy + CAI
    // cost. Bounded by 4^size paths. Keeps the DP cost consistent with the structure evaluator.
    void hp_enum(int cur, int bp, int size, int nt5, int nt3, std::string& loop, double cai,
                 double& best, std::string* best_loop) {
        if ((int)loop.size() == size) {
            if (cur != bp) return;
            std::string span; span.reserve(size + 2);
            span.push_back(nuc_char(nt5)); span += loop; span.push_back(nuc_char(nt3));
            int si = energy::special_index(span);
            int e = energy::hairpin(size, nt5, nuc_code(loop.front()), nuc_code(loop.back()), nt3, si);
            double v = (double)e + cai;
            if (v < best) { best = v; if (best_loop) *best_loop = loop; }
            return;
        }
        for (const Edge& e : d_.out_edges(cur)) {
            loop.push_back(nuc_char(e.nuc));
            hp_enum(e.to, bp, size, nt5, nt3, loop, cai + lambda_ * e.weight, best, best_loop);
            loop.pop_back();
        }
    }

    double gap_min(int x, int y) {
        double best = kInf;
        for (auto& g : gap_options(x, y)) best = std::min(best, g.cai);
        return best;
    }
    // reconstruct the nucleotides of an unpaired stretch matching a chosen gap option
    std::string gap_string(int x, int y, const GapOpt& g) {
        int dp = d_.pos_of(y) - d_.pos_of(x);
        if (dp == 0) return "";
        if (dp == 1) return std::string(1, nuc_char(g.first));
        return nuc_char(g.first) + mincai_path(g.x2, g.y2) + nuc_char(g.last);
    }
    std::string mincai_path(int a, int b) {
        std::string s;
        while (a != b) {
            double tot = mincai(a, b);
            bool done = false;
            for (const Edge& e : d_.out_edges(a)) {
                double sub = mincai(e.to, b);
                if (sub < kInf && eq(lambda_ * e.weight + sub, tot)) {
                    s.push_back(nuc_char(e.nuc)); a = e.to; done = true; break;
                }
            }
            if (!done) break;
        }
        return s;
    }

    // two-loop (stack/bulge/internal) closed by (nt5@i, nt3@j-1), inner span ap..bp.
    double two_loop_cost(int ap, int bp, int nt5, int nt3, TLDetail* det) {
        double best = kInf;
        int posap = d_.pos_of(ap), posbp = d_.pos_of(bp);
        // The outer closing pair (nt5,nt3) never changes within this call, so its pair-type is a
        // whole-function invariant; energy::two_loop() used to recompute pair_type(nt5,nt3) on
        // every single (p,ia,q1,ib,q[,a5,a3]) combination below (millions of calls for a large
        // protein) via the un-hoisted two_loop() wrapper. Compute it once and call two_loop_pt().
        int type_outer = energy::pair_type(nt5, nt3);
        for (int p = posap; p <= posbp - 2; ++p) {
            int n1 = p - posap;
            if (n1 > kMaxLoop) break;
            for (int ia : d_.nodes_at(p)) {
                if (n1 == 0 && ia != ap) continue;
                // inner pair span >= min_loop+2; total unpaired n1+n2 <= kMaxLoop.
                int q1lo = std::max(p + min_loop_ + 2, posbp - (kMaxLoop - n1));
                for (int q1 = q1lo; q1 <= posbp; ++q1) {
                    int n2 = posbp - q1;
                    if (n2 < 0 || n1 + n2 > kMaxLoop) continue;
                    for (int ib : d_.nodes_at(q1)) {
                        if (n2 == 0 && ib != bp) continue;
                        const auto& inner = closed(ia, ib);
                        if (inner.empty()) continue;
                        if (n1 == 0 && n2 == 0) {                       // stack
                            for (const PairCost& q : inner) {
                                int type2 = energy::pair_type(q.nt3, q.nt5);
                                int e = energy::two_loop_pt(0, 0, type_outer, nt5, nt3, nt5, type2, nt3);
                                double v = (double)e + q.cost;
                                if (v < best) { best = v;
                                    if (det) *det = TLDetail{ia, ib, q.nt5, q.nt3, -1, -1, -1, -1, {}, {}}; }
                            }
                        } else if (n1 == 0 || n2 == 0) {                // bulge
                            const GapOpt* bg = bulge_gap(n1 == 0 ? ib : ap, n1 == 0 ? bp : ia);
                            if (!bg) continue;
                            for (const PairCost& q : inner) {
                                int type2 = energy::pair_type(q.nt3, q.nt5);
                                int e = energy::two_loop_pt(n1, n2, type_outer, nt5, nt3, nt5, type2, nt3);
                                double v = (double)e + q.cost + bg->cai;
                                if (v < best) { best = v;
                                    if (det) {
                                        // bulge unpaired is on the 3' side when n1==0, else the 5' side.
                                        if (n1 == 0) *det = TLDetail{ia, ib, q.nt5, q.nt3, -1, -1, ib, bp, {}, *bg};
                                        else         *det = TLDetail{ia, ib, q.nt5, q.nt3, ap, ia, -1, -1, *bg, {}};
                                    } }
                            }
                        } else {                                        // internal loop
                            const auto& g5 = gap_options(ap, ia);
                            const auto& g3 = gap_options(ib, bp);
                            for (const PairCost& q : inner) {
                                // type2 depends only on q, not on a5/a3: hoist out of their double loop.
                                int type2 = energy::pair_type(q.nt3, q.nt5);
                                for (auto& a5 : g5)
                                    for (auto& a3 : g3) {
                                        int e = energy::two_loop_pt(n1, n2, type_outer, a5.first, a3.last,
                                                                    a5.last, type2, a3.first);
                                        double v = (double)e + q.cost + a5.cai + a3.cai;
                                        if (v < best) { best = v;
                                            if (det) *det = TLDetail{ia, ib, q.nt5, q.nt3,
                                                                     ap, ia, ib, bp, a5, a3}; }
                                    }
                            }
                        }
                    }
                }
            }
        }
        return best;
    }
    const GapOpt* bulge_gap(int x, int y) {
        const auto& gs = gap_options(x, y);
        const GapOpt* best = nullptr;
        for (auto& g : gs) if (!best || g.cai < best->cai) best = &g;
        return best;
    }

    double m1(int a, int b) {
        if (a == b) return kInf;
        size_t k = key(a, b);
        if (m1_v_[k] != kUncomputed) return m1_v_[k];
        double best = kInf;
        for (const Edge& e : d_.out_edges(a)) {
            double sub = m1(e.to, b);
            if (sub < kInf) { double v = lambda_ * e.weight + energy::multi_unpaired() + sub;
                if (v < best) best = v; }
        }
        each_branch(a, b, [&](double branch, int c, int n5, int n3) {
            double restu = multi_unpaired_gap(c, b);
            if (restu < kInf) {
                double v = branch + restu;
                if (v < best) best = v; }
            double more = m1(c, b);
            if (more < kInf) { double v = branch + more;
                if (v < best) best = v; }
        });
        m1_v_[k] = best;
        return best;
    }
    double m2(int a, int b) {
        if (a == b) return kInf;
        size_t k = key(a, b);
        if (m2_v_[k] != kUncomputed) return m2_v_[k];
        double best = kInf;
        for (const Edge& e : d_.out_edges(a)) {
            double sub = m2(e.to, b);
            if (sub < kInf) { double v = lambda_ * e.weight + energy::multi_unpaired() + sub;
                if (v < best) best = v; }
        }
        each_branch(a, b, [&](double branch, int c, int n5, int n3) {
            double more = m1(c, b);
            if (more < kInf) { double v = branch + more;
                if (v < best) best = v; }
        });
        m2_v_[k] = best;
        return best;
    }

    // Fused dense m1/m2 kernel for the validation wavefront. Both states enumerate the same eager
    // branch row, so one pass preserves each state's original candidate order while avoiding a
    // duplicate scan. At the current gap brow_[a] contains exactly the branches ending no later
    // than b's nucleotide position, in each_branch() order.
    template <bool Profile>
    void m1m2_dense(int a, int b) {
        const size_t k = key(a, b);
        double best1 = kInf;
        double best2 = kInf;
        if (a != b) {
            const int pa = d_.pos_of(a);
            const int pb = d_.pos_of(b);
            // A one-nucleotide interval has no non-empty remainder after consuming its first edge.
            if (pb - pa > 1) {
                for (const Edge& e : d_.out_edges(a)) {
                    double sub1 = m1(e.to, b);
                    if (sub1 < kInf) {
                        double v = lambda_ * e.weight + energy::multi_unpaired() + sub1;
                        if (v < best1) best1 = v;
                    }
                    double sub2 = m2(e.to, b);
                    if (sub2 < kInf) {
                        double v = lambda_ * e.weight + energy::multi_unpaired() + sub2;
                        if (v < best2) best2 = v;
                    }
                }
            }
            for (const BranchEntry& branch : brow_[a]) {
                if (Profile) ++split_visits_[a];
                const int pc = d_.pos_of(branch.c);
                if (pc > pb) break;
                if (branch.c == b) {
                    // Empty all-unpaired suffix after one direct branch.
                    if (branch.bcost < best1) best1 = branch.bcost;
                    continue;
                }
                // Distinct DFA states at the same layer have no epsilon path between them.
                if (pc == pb) continue;
                if (Profile) ++split_eligible_[a];
                double restu = multi_unpaired_gap(branch.c, b);
                if (restu < kInf) {
                    double v = branch.bcost + restu;
                    if (v < best1) best1 = v;
                }
                double more = m1(branch.c, b);
                if (more < kInf) {
                    double v = branch.bcost + more;
                    if (v < best1) best1 = v;
                    if (v < best2) best2 = v;
                }
            }
        }
        m1_v_[k] = best1;
        m2_v_[k] = best2;
    }

    // Exact sparse Turner multiloop recurrence. A multiloop fragment is
    // either extended by an unpaired edge, is one direct closed branch, or is partitionable into
    // at least two branches. For m2(a,b), it is sufficient to consider rightmost direct branches
    // [c,b] that beat every partitionable/unpaired alternative for their own interval. If a direct
    // branch is not such a candidate, reassociating through its better alternative is no worse
    // (the standard sparse-MFE triangle-inequality argument). Worst-case time is still cubic, but
    // practical time is O(n^2 + nZ), where Z is the number of retained candidates.
    template <bool Prune, bool Profile>
    void m1m2_right_normal(int a, int b) {
        size_t k = key(a, b);
        double best2 = kInf;
        if (a != b) {
            // Unpaired suffix after a fragment with at least two branches.
            for (const Edge& e : d_.in_edges(b)) {
                double sub = m2(a, e.from);
                if (sub < kInf) {
                    double v = sub + lambda_ * e.weight + energy::multi_unpaired();
                    if (v < best2) best2 = v;
                }
            }
            // Partition at the start of a sparse rightmost branch candidate [c,b]. Candidates
            // at the same nucleotide position as a cannot follow a non-empty left fragment.
            int pa = d_.pos_of(a);
            for (const MultiCandidate& cand : mcand_end_[b]) {
                if (Profile) ++split_visits_[b];
                if (d_.pos_of(cand.start) <= pa) continue;
                if (Profile) ++split_eligible_[b];
                double left = m1(a, cand.start);
                if (left < kInf) {
                    double v = left + cand.bcost;
                    if (v < best2) best2 = v;
                }
            }
        }
        m2_v_[k] = best2;

        // Best partitionable/unpaired one-or-more-branch fragment. Both endpoint extensions are
        // standard even for a fixed sequence. On a codon DFA we must additionally range over every
        // actual incident edge and keep right-branch candidates keyed by the concrete endpoint
        // state, because different paths can branch or merge at either side.
        double wp = best2;
        if (a != b) {
            for (const Edge& e : d_.out_edges(a)) {
                double sub = m1(e.to, b);
                if (sub < kInf) {
                    double v = lambda_ * e.weight + energy::multi_unpaired() + sub;
                    if (v < wp) wp = v;
                }
            }
            for (const Edge& e : d_.in_edges(b)) {
                double sub = m1(a, e.from);
                if (sub < kInf) {
                    double v = sub + lambda_ * e.weight + energy::multi_unpaired();
                    if (v < wp) wp = v;
                }
            }
        }

        // Direct branch V(a,b)+b: collapse pair types because every downstream recurrence sees
        // only its scalar cost. Strict '<' retains the original first-pair tie order.
        MultiCandidate direct{a, kInf};
        for (const PairCost& pc : closed(a, b)) {
            double v = pc.cost + energy::multi_branch(pc.nt5, pc.nt3);
            if (v < direct.bcost)
                direct = MultiCandidate{a, v};
        }
        m1_v_[k] = std::min(wp, direct.bcost);

        // Candidate criterion. A tie is partitionable and need not be retained; backtrace still
        // replays the original dense recurrence from the completed scalar costs.
        if (direct.bcost < (Prune ? wp : kInf)) mcand_end_[b].push_back(direct);
    }
    // Backtrace-only replay of the original left-to-right branch enumeration. Rows are extended
    // lazily and remain position-sorted, so reconstruction stores only the small subset it visits.
    template <class F>
    void each_branch(int a, int b, F f) {
        int pa = d_.pos_of(a), pb = d_.pos_of(b);
        int& hi = brow_hi_[a];
        if (pb > hi) {
            int start = (hi < 0) ? pa + min_loop_ + 2 : hi + 1;
            for (int m = start; m <= pb; ++m)
                for (int c : d_.nodes_at(m))
                    for (const PairCost& pc : closed(a, c))
                        brow_[a].push_back(BranchEntry{c, pc.nt5, pc.nt3,
                                            pc.cost + energy::multi_branch(pc.nt5, pc.nt3)});
            hi = pb;
        }
        for (const BranchEntry& e : brow_[a]) {
            if (d_.pos_of(e.c) > pb) break;
            f(e.bcost, e.c, e.nt5, e.nt3);
        }
    }
    const std::vector<PairCost>& closed(int a, int b) {
        size_t k = key(a, b);
        if (closed_c_[k]) return closed_v_[k];
        // Build into a local vector and commit at the end (mirrors the original): closed(a,b)
        // only ever recurses into strictly-smaller spans (the DFA is position-increasing), so
        // it never re-enters its own cell mid-computation.
        std::vector<PairCost> res;
        if (d_.pos_of(b) - d_.pos_of(a) >= min_loop_ + 2) {
            // Several (out-edge, in-edge) pairs can share the same closing pair (nt5,nt3) but
            // reach different inner spans (the DFA's tier structure). Keep the argmin per
            // (nt5,nt3).
            // There are only six pairable nucleotide pairs. A fixed table avoids constructing a
            // red-black tree (and allocating up to six nodes) for every closed cell. Iterating the
            // numeric keys below preserves std::map's original result/tie order exactly.
            double bestpair[40];
            std::fill(bestpair, bestpair + 40, kInf);
            for (const Edge& ea : d_.out_edges(a))
                for (const Edge& eb : d_.in_edges(b)) {
                    if (!energy::pairable(ea.nuc, eb.nuc)) continue;
                    int ap = ea.to, bp = eb.from;
                    double hc = hairpin_cost(ap, bp, ea.nuc, eb.nuc, nullptr);
                    double tc = two_loop_cost(ap, bp, ea.nuc, eb.nuc, nullptr);
                    double mc = m2(ap, bp);
                    double mco = (mc < kInf) ? energy::multi_closing(ea.nuc, eb.nuc) + mc : kInf;
                    double loop = std::min(hc, std::min(tc, mco));
                    if (loop >= kInf) continue;
                    double cost = lambda_ * (ea.weight + eb.weight) + loop;
                    int pk = ea.nuc * 8 + eb.nuc;
                    if (cost < bestpair[pk]) bestpair[pk] = cost;
                }
            int pair_count = 0;
            for (int pk = 0; pk < 40; ++pk)
                if (bestpair[pk] < kInf) ++pair_count;
            res.reserve(pair_count);
            for (int pk = 0; pk < 40; ++pk) {
                if (bestpair[pk] >= kInf) continue;
                res.push_back(PairCost{pk / 8, pk % 8, bestpair[pk]});
            }
        }
        closed_v_[k] = std::move(res); closed_c_[k] = 1;
        return closed_v_[k];
    }

    double ext(int a, int b) {
        if (a == b) return 0.0;
        if (ext_v_[a] != kUncomputed) return ext_v_[a];
        double best = kInf;
        for (const Edge& e : d_.out_edges(a)) {
            double sub = ext(e.to, b);
            if (sub < kInf) { double v = lambda_ * e.weight + energy::external_unpaired() + sub;
                if (v < best) best = v; }
        }
        int pa = d_.pos_of(a), pb = d_.pos_of(b);
        for (int m = pa + min_loop_ + 2; m <= pb; ++m)
            for (int c : d_.nodes_at(m))
                for (const PairCost& pc : closed(a, c)) {
                    double sub = ext(c, b);
                    if (sub < kInf) { double v = pc.cost + energy::external_paired(pc.nt5, pc.nt3) + sub;
                        if (v < best) best = v; }
                }
        ext_v_[a] = best;
        return best;
    }

    // ---------------- backtrace ----------------
    // Decisions are intentionally reconstructed from the completed cost memos instead of stored
    // for every O(n^2) cell. Only O(n) decisions on the winning derivation are needed, and replaying
    // the original candidate order preserves the exact tie-break while avoiding large decision-table
    // allocations and write traffic during the hot forward pass.
    NodeDec ext_decision(int a, int b) {
        double best = kInf;
        NodeDec dec{-1, 0, 0, 0, 0, 0};
        for (const Edge& e : d_.out_edges(a)) {
            double sub = ext(e.to, b);
            if (sub < kInf) {
                double v = lambda_ * e.weight + energy::external_unpaired() + sub;
                if (v < best) { best = v; dec = NodeDec{0, e.nuc, e.to, 0, 0, 0}; }
            }
        }
        int pa = d_.pos_of(a), pb = d_.pos_of(b);
        for (int m = pa + min_loop_ + 2; m <= pb; ++m)
            for (int c : d_.nodes_at(m))
                for (const PairCost& pc : closed(a, c)) {
                    double sub = ext(c, b);
                    if (sub < kInf) {
                        double v = pc.cost + energy::external_paired(pc.nt5, pc.nt3) + sub;
                        if (v < best) {
                            best = v;
                            dec = NodeDec{1, 0, 0, c, pc.nt5, pc.nt3};
                        }
                    }
                }
        return dec;
    }

    NodeDec m1_decision(int a, int b) {
        double best = kInf;
        NodeDec dec{-1, 0, 0, 0, 0, 0};
        for (const Edge& e : d_.out_edges(a)) {
            double sub = m1(e.to, b);
            if (sub < kInf) {
                double v = lambda_ * e.weight + energy::multi_unpaired() + sub;
                if (v < best) { best = v; dec = NodeDec{0, e.nuc, e.to, 0, 0, 0}; }
            }
        }
        each_branch(a, b, [&](double branch, int c, int nt5, int nt3) {
            double restu = multi_unpaired_gap(c, b);
            if (restu < kInf) {
                double v = branch + restu;
                if (v < best) { best = v; dec = NodeDec{2, 0, 0, c, nt5, nt3}; }
            }
            double more = m1(c, b);
            if (more < kInf) {
                double v = branch + more;
                if (v < best) { best = v; dec = NodeDec{3, 0, 0, c, nt5, nt3}; }
            }
        });
        return dec;
    }

    NodeDec m2_decision(int a, int b) {
        double best = kInf;
        NodeDec dec{-1, 0, 0, 0, 0, 0};
        for (const Edge& e : d_.out_edges(a)) {
            double sub = m2(e.to, b);
            if (sub < kInf) {
                double v = lambda_ * e.weight + energy::multi_unpaired() + sub;
                if (v < best) { best = v; dec = NodeDec{0, e.nuc, e.to, 0, 0, 0}; }
            }
        }
        each_branch(a, b, [&](double branch, int c, int nt5, int nt3) {
            double more = m1(c, b);
            if (more < kInf) {
                double v = branch + more;
                if (v < best) { best = v; dec = NodeDec{3, 0, 0, c, nt5, nt3}; }
            }
        });
        return dec;
    }

    ClosedDec closed_decision(int a, int b, int nt5, int nt3) {
        double best = kInf;
        ClosedDec dec{};
        for (const Edge& ea : d_.out_edges(a))
            for (const Edge& eb : d_.in_edges(b)) {
                if (ea.nuc != nt5 || eb.nuc != nt3) continue;
                int ap = ea.to, bp = eb.from;
                HpDec hp;
                TLDetail tl;
                double hc = hairpin_cost(ap, bp, nt5, nt3, &hp);
                double tc = two_loop_cost(ap, bp, nt5, nt3, &tl);
                double mc = m2(ap, bp);
                double mco = (mc < kInf) ? energy::multi_closing(nt5, nt3) + mc : kInf;
                double loop = std::min(hc, std::min(tc, mco));
                if (loop >= kInf) continue;
                double cost = lambda_ * (ea.weight + eb.weight) + loop;
                if (cost >= best) continue;
                best = cost;
                dec = ClosedDec{ap, bp, 0, std::move(hp), tl};
                if (eq(loop, hc)) dec.kind = 0;
                else if (eq(loop, tc)) dec.kind = 1;
                else dec.kind = 2;
            }
        return dec;
    }

    void bt_ext(int a, int b, std::string& seq, std::string& st) {
        if (a == b) return;
        NodeDec d = ext_decision(a, b);
        if (d.kind == 0) { seq.push_back(nuc_char(d.nt)); st.push_back('.'); bt_ext(d.to, b, seq, st); }
        else { bt_closed(a, d.c, d.nt5, d.nt3, seq, st); bt_ext(d.c, b, seq, st); }
    }
    void bt_m1(int a, int b, std::string& seq, std::string& st) {
        NodeDec d = m1_decision(a, b);
        if (d.kind == 0) { seq.push_back(nuc_char(d.nt)); st.push_back('.'); bt_m1(d.to, b, seq, st); }
        else if (d.kind == 2) { bt_closed(a, d.c, d.nt5, d.nt3, seq, st);
            seq += mincai_path(d.c, b); st.append(d_.pos_of(b) - d_.pos_of(d.c), '.'); }
        else { bt_closed(a, d.c, d.nt5, d.nt3, seq, st); bt_m1(d.c, b, seq, st); }
    }
    void bt_m2(int a, int b, std::string& seq, std::string& st) {
        NodeDec d = m2_decision(a, b);
        if (d.kind == 0) { seq.push_back(nuc_char(d.nt)); st.push_back('.'); bt_m2(d.to, b, seq, st); }
        else { bt_closed(a, d.c, d.nt5, d.nt3, seq, st); bt_m1(d.c, b, seq, st); }
    }
    void bt_closed(int a, int b, int nt5, int nt3, std::string& seq, std::string& st) {
        ClosedDec d = closed_decision(a, b, nt5, nt3);
        seq.push_back(nuc_char(nt5)); st.push_back('(');
        if (d.kind == 0) {                     // hairpin
            std::string loop = d.hp.special ? d.hp.loop : gap_string(d.ap, d.bp, d.hp.g);
            seq += loop; st.append(loop.size(), '.');
        } else if (d.kind == 1) {              // two-loop
            std::string g5 = (d.tl.x5 >= 0) ? gap_string(d.tl.x5, d.tl.y5, d.tl.g5) : "";
            std::string g3 = (d.tl.x3 >= 0) ? gap_string(d.tl.x3, d.tl.y3, d.tl.g3) : "";
            seq += g5; st.append(g5.size(), '.');
            bt_closed(d.tl.ia, d.tl.ib, d.tl.ntp, d.tl.ntq, seq, st);
            seq += g3; st.append(g3.size(), '.');
        } else {                               // multiloop
            bt_m2(d.ap, d.bp, seq, st);
        }
        seq.push_back(nuc_char(nt3)); st.push_back(')');
    }

    const DFA& d_;
    double lambda_;
    int min_loop_;
    static void record_phase(double& seconds, std::chrono::steady_clock::time_point& start) {
        const auto now = std::chrono::steady_clock::now();
        seconds += std::chrono::duration<double>(now - start).count();
        start = now;
    }
    TurnerProfile profile_;
    // Updated by the unique owner of an endpoint/start row in each wavefront phase.
    std::vector<size_t> split_visits_, split_eligible_;

    int n_ = 0;  // num_nodes(), cached for flat indexing
    size_t memo_cells_ = 0;
#ifndef LDCLEAN_SQUARE_MEMOS
    std::vector<size_t> row_base_;       // adjusted packed-row offset, indexed by node ID
#endif
#ifdef LDCLEAN_PACKED_RANK_MAP
    std::vector<int> rank_by_position_;  // validation C: redundant node ID -> packed rank map
#endif
    int num_threads_ = 1;

    // Dense scalar memo tables. The production layout is square; the low-memory build stores only
    // layer-forward cells. kUncomputed is below all valid costs; +infinity means infeasible.
    std::vector<double> mincai_v_, m1_v_, m2_v_;
    std::vector<double> ext_v_;  // 1-D: ext's second node is always d_.final().

    // Dense vector-valued memos keyed by (a,b): the closed-pair list and the unpaired-gap
    // option list. These are the two hottest tables after the scalar memos (closed() alone
    // accounts for the majority of calls), so densifying them removes their hash overhead too.
    std::vector<std::vector<PairCost>> closed_v_;
    std::vector<std::vector<GapOpt>>   gap_v_;
    std::vector<char>                  closed_c_, gap_c_;

    // Dense branch rows materialized only for states visited by backtrace. brow_hi_ records the
    // highest covered position, allowing a later query to extend rather than rebuild a row.
    std::vector<std::vector<BranchEntry>> brow_;
    std::vector<int> brow_hi_;
    // Exact sparse candidates keyed by right endpoint, used only by the bottom-up path. Each
    // candidate replaces all non-winning pair types for its cell in the multiloop bifurcation.
    std::vector<std::vector<MultiCandidate>> mcand_end_;

};

}  // namespace

int turner_effective_threads(int requested) {
#ifdef _OPENMP
    const int capped = std::max(1, std::min(requested, useful_thread_ceiling()));
    // `num_threads(capped)` is still only a request: a non-OpenMP build, OMP_THREAD_LIMIT,
    // OMP_DYNAMIC, or a surrounding parallel region may supply a smaller team. Probe the same
    // request used by the solver so callers do not report the policy cap as an observed team size.
    int observed = 1;
#pragma omp parallel num_threads(capped)
    {
#pragma omp single
        observed = omp_get_num_threads();
    }
    return observed;
#else
    (void)requested;
    return 1;
#endif
}

double turner_min_cost(const DFA& dfa, double lambda, int min_loop, int threads,
                       TurnerSparseStats* sparse_stats) {
    TurnerSolver solver(dfa, lambda, min_loop, threads);
    double cost = solver.solve();
    if (sparse_stats) solver.sparse_stats(*sparse_stats);
    return cost;
}

double turner_min_cost_dense_reference(const DFA& dfa, double lambda, int min_loop) {
    return TurnerSolver(dfa, lambda, min_loop, 1).solve_dense_reference();
}

double turner_min_cost_dense_wavefront(const DFA& dfa, double lambda, int min_loop, int threads) {
    return TurnerSolver(dfa, lambda, min_loop, threads).solve_dense_wavefront();
}

double turner_min_cost_dense_right_normal(const DFA& dfa, double lambda, int min_loop,
                                          int threads, TurnerSparseStats* sparse_stats) {
    TurnerSolver solver(dfa, lambda, min_loop, threads);
    double cost = solver.solve_dense_right_normal();
    if (sparse_stats) solver.sparse_stats(*sparse_stats);
    return cost;
}

double turner_min_cost_profiled(const DFA& dfa, double lambda, TurnerKernel kernel,
                                TurnerProfile& profile, int min_loop, int threads,
                                TurnerSparseStats* sparse_stats) {
    profile = TurnerProfile{};
    const auto start = std::chrono::steady_clock::now();
    TurnerSolver solver(dfa, lambda, min_loop, threads);
    profile.construction_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    double cost = solver.solve_profiled(kernel, profile);
    if (sparse_stats) solver.sparse_stats(*sparse_stats);
    return cost;
}

TurnerDesign design_turner(const DFA& dfa, const CodonTable& table, double lambda, int min_loop,
                           int threads, TurnerSparseStats* sparse_stats) {
    TurnerSolver solver(dfa, lambda, min_loop, threads);
    TurnerDesign r;
    r.cost = solver.solve();
    if (!std::isfinite(r.cost))
        throw std::overflow_error("Turner objective exceeds the finite double range");
    if (sparse_stats) solver.sparse_stats(*sparse_stats);
    solver.backtrace(r.mrna, r.structure);
    r.mfe_kcal = turner_eval(r.mrna, r.structure) / 100.0;
    r.cai = table.cai(r.mrna);
    r.codon_penalty = 0.0;
    for (size_t i = 0; i < r.mrna.size(); i += 3)
        r.codon_penalty -= std::log(table.w_of(r.mrna.substr(i, 3)));
    // Independently recover the selected path's edge weights. This also correctly handles the
    // validation constructor's zero-weight fixed RNA, for which the table penalty is unrelated.
    int node = dfa.start();
    double path_cost = 0.0;
    for (char nucleotide : r.mrna) {
        const Edge* chosen = nullptr;
        for (const Edge& edge : dfa.out_edges(node))
            if (edge.nuc == nuc_code(nucleotide)) {
                if (chosen) throw std::logic_error("ambiguous reconstructed DFA path");
                chosen = &edge;
            }
        if (!chosen) throw std::logic_error("reconstructed RNA leaves the DFA");
        path_cost += chosen->weight;
        node = chosen->to;
    }
    if (node != dfa.final()) throw std::logic_error("reconstructed RNA does not reach DFA final state");
    r.realized_cost = r.mfe_kcal * 100.0 + lambda * path_cost;
    const double tolerance = 1e-4 + 1e-10 * std::max(std::fabs(r.cost), std::fabs(r.realized_cost));
    if (!std::isfinite(r.realized_cost) || std::fabs(r.cost - r.realized_cost) > tolerance)
        throw std::logic_error("reconstructed objective differs from the dynamic-programming optimum");
    return r;
}

// ---------------- independent structure-energy evaluator (Turner d0) ----------------
namespace {
// energy of the loop closed by pair (i,j), recursing into nested children.
int eval_closed(const std::string& seq, const std::vector<int>& pair, int i, int j) {
    int nuci = nuc_code(seq[i]), nucj = nuc_code(seq[j]);
    // collect immediate children pairs
    std::vector<std::pair<int, int>> kids;
    for (int k = i + 1; k < j;) {
        if (pair[k] > k) { kids.emplace_back(k, pair[k]); k = pair[k] + 1; }
        else ++k;
    }
    if (kids.empty()) {  // hairpin
        int size = j - i - 1;
        std::string span = seq.substr(i, j - i + 1);
        int si = energy::special_index(span);
        return energy::hairpin(size, nuci, nuc_code(seq[i + 1]), nuc_code(seq[j - 1]), nucj, si);
    }
    if (kids.size() == 1) {  // two-loop
        int p = kids[0].first, q = kids[0].second;
        int e = energy::two_loop(p - i - 1, j - q - 1, nuci, nuc_code(seq[i + 1]),
                                 nuc_code(seq[j - 1]), nucj, nuc_code(seq[p - 1]), nuc_code(seq[p]),
                                 nuc_code(seq[q]), nuc_code(seq[q + 1]));
        return e + eval_closed(seq, pair, p, q);
    }
    // multiloop
    int e = energy::multi_closing(nuci, nucj);
    int unpaired = (j - i - 1);
    for (auto& kv : kids) {
        e += energy::multi_branch(nuc_code(seq[kv.first]), nuc_code(seq[kv.second]));
        e += eval_closed(seq, pair, kv.first, kv.second);
        unpaired -= (kv.second - kv.first + 1);
    }
    e += unpaired * energy::multi_unpaired();
    return e;
}
}  // namespace

int turner_eval(const std::string& seq, const std::string& structure) {
    energy::init();
    int n = (int)seq.size();
    std::vector<int> pair(n, -1), stack;
    for (int k = 0; k < n; ++k) {
        if (structure[k] == '(') stack.push_back(k);
        else if (structure[k] == ')') { int o = stack.back(); stack.pop_back(); pair[o] = k; pair[k] = o; }
    }
    // external loop: sum external stems (+ external_paired) and unpaired (0)
    int e = 0;
    for (int k = 0; k < n;) {
        if (pair[k] > k) {
            e += energy::external_paired(nuc_code(seq[k]), nuc_code(seq[pair[k]]));
            e += eval_closed(seq, pair, k, pair[k]);
            k = pair[k] + 1;
        } else { e += energy::external_unpaired(); ++k; }
    }
    return e;
}

}  // namespace ldclean
