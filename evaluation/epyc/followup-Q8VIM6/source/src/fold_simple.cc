// fold_simple.cc — see fold_simple.h.
#include "fold_simple.h"

#include <limits>
#include <unordered_map>
#include <vector>

namespace ldclean {

bool simple_pairable(int a, int b) {
    return (a == 2 && b == 3) || (a == 3 && b == 2) ||  // C-G / G-C
           (a == 1 && b == 4) || (a == 4 && b == 1) ||  // A-U / U-A
           (a == 3 && b == 4) || (a == 4 && b == 3);    // G-U / U-G
}
int simple_pair_score(int a, int b) {
    if ((a == 2 && b == 3) || (a == 3 && b == 2)) return -3;  // C-G
    if ((a == 1 && b == 4) || (a == 4 && b == 1)) return -2;  // A-U
    if ((a == 3 && b == 4) || (a == 4 && b == 3)) return -1;  // G-U
    return 0;
}

namespace {

constexpr double kInf = 1e18;

// Plain aggregates (no default member initializers, to keep C++11 brace-init working).
struct SBack { int type; int next; int nt; };  // type 0: unpaired (next=a', nt); 1: bifurc (next=c)
struct PBack { int ap; int bp; int nta; int ntb; };

class SimpleSolver {
public:
    SimpleSolver(const DFA& d, double lambda, int min_loop)
        : d_(d), lambda_(lambda), min_span_(min_loop + 2) {}

    double solve_S(int a, int b) {
        if (a == b) return 0.0;
        long long key = pack(a, b);
        auto it = s_.find(key);
        if (it != s_.end()) return it->second;
        s_[key] = kInf;  // guard against re-entry (DAG, so this is just memo init)

        double best = kInf;
        SBack bk{};
        // first nucleotide unpaired
        for (const Edge& e : d_.out_edges(a)) {
            double sub = solve_S(e.to, b);
            if (sub >= kInf) continue;
            double v = lambda_ * e.weight + sub;
            if (v < best) { best = v; bk = SBack{0, e.to, e.nuc}; }
        }
        // first block is a pair (bifurcation S -> P S)
        int pa = d_.pos_of(a), pb = d_.pos_of(b);
        for (int m = pa + min_span_; m <= pb; ++m) {
            for (int c : d_.nodes_at(m)) {
                double pp = solve_P(a, c);
                if (pp >= kInf) continue;
                double ss = solve_S(c, b);
                if (ss >= kInf) continue;
                double v = pp + ss;
                if (v < best) { best = v; bk = SBack{1, c, 0}; }
            }
        }
        s_[key] = best;
        sback_[key] = bk;
        return best;
    }

    double solve_P(int a, int b) {
        if (d_.pos_of(b) - d_.pos_of(a) < min_span_) return kInf;
        long long key = pack(a, b);
        auto it = p_.find(key);
        if (it != p_.end()) return it->second;
        p_[key] = kInf;

        double best = kInf;
        PBack bk{};
        for (const Edge& ea : d_.out_edges(a)) {          // a -(nta)-> a'  (5' nt at pos i)
            for (const Edge& eb : d_.in_edges(b)) {        // b' -(ntb)-> b  (3' nt at pos j-1)
                if (!simple_pairable(ea.nuc, eb.nuc)) continue;
                double inner = solve_S(ea.to, eb.from);
                if (inner >= kInf) continue;
                double v = lambda_ * (ea.weight + eb.weight) +
                           simple_pair_score(ea.nuc, eb.nuc) + inner;
                if (v < best) { best = v; bk = PBack{ea.to, eb.from, ea.nuc, eb.nuc}; }
            }
        }
        p_[key] = best;
        pback_[key] = bk;
        return best;
    }

    void back_S(int a, int b, std::string& seq, std::string& st) {
        if (a == b) return;
        const SBack& bk = sback_[pack(a, b)];
        if (bk.type == 0) {
            seq.push_back(nuc_char(bk.nt));
            st.push_back('.');
            back_S(bk.next, b, seq, st);
        } else {
            back_P(a, bk.next, seq, st);
            back_S(bk.next, b, seq, st);
        }
    }
    void back_P(int a, int b, std::string& seq, std::string& st) {
        const PBack& bk = pback_[pack(a, b)];
        seq.push_back(nuc_char(bk.nta));
        st.push_back('(');
        back_S(bk.ap, bk.bp, seq, st);
        seq.push_back(nuc_char(bk.ntb));
        st.push_back(')');
    }

private:
    long long pack(int a, int b) const { return (long long)a * d_.num_nodes() + b; }

    const DFA& d_;
    double lambda_;
    int min_span_;
    std::unordered_map<long long, double> s_, p_;
    std::unordered_map<long long, SBack> sback_;
    std::unordered_map<long long, PBack> pback_;
};

}  // namespace

DesignResult design_simple(const DFA& dfa, double lambda, int min_loop) {
    SimpleSolver solver(dfa, lambda, min_loop);
    double cost = solver.solve_S(dfa.start(), dfa.final());

    DesignResult r;
    r.cost = cost;
    solver.back_S(dfa.start(), dfa.final(), r.mrna, r.structure);

    // Split the cost into folding energy (sum of pair scores) and CAI cost.
    double fold = 0.0;
    std::vector<int> stack;
    for (size_t k = 0; k < r.structure.size(); ++k) {
        if (r.structure[k] == '(') stack.push_back((int)k);
        else if (r.structure[k] == ')') {
            int i = stack.back(); stack.pop_back();
            fold += simple_pair_score(nuc_code(r.mrna[i]), nuc_code(r.mrna[k]));
        }
    }
    r.fold_energy = fold;
    r.cai_cost = cost - fold;
    return r;
}

}  // namespace ldclean
