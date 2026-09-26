// Exact-integer, all-cell validation of the abstract multiloop sparsification kernel.
//
// This test is intentionally independent of RNA energy code. It generates small layered DAGs with
// branching, merging, unreachable state pairs, negative edge costs, and arbitrary direct-branch
// hyperedges. It compares, for every interval:
//   * exhaustive shortest hyperpaths with >=1 and >=2 branch atoms;
//   * a dense left-normal recurrence (the shape retained as the C++ validation oracle);
//   * a dense right-normal recurrence; and
//   * the candidate-restricted right-normal recurrence used by the production solver.
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

namespace {

using Cost = std::int64_t;
constexpr Cost INF = (Cost{1} << 60);

Cost add(Cost a, Cost b) {
    return (a >= INF || b >= INF) ? INF : a + b;
}

struct Arc { int from, to; Cost cost; };

struct Graph {
    int layers = 0;
    std::vector<int> pos;
    std::vector<std::vector<int>> at;
    std::vector<std::vector<Arc>> out, in;
    std::vector<Cost> branch;  // flattened B(a,b)

    int nodes() const { return static_cast<int>(pos.size()); }
    Cost B(int a, int b) const { return branch[static_cast<size_t>(a) * nodes() + b]; }
};

struct Tables {
    int n = 0;
    std::vector<Cost> m1, m2;
    explicit Tables(int nodes) : n(nodes), m1(static_cast<size_t>(n) * n, INF),
                                 m2(static_cast<size_t>(n) * n, INF) {}
    Cost& M1(int a, int b) { return m1[static_cast<size_t>(a) * n + b]; }
    Cost& M2(int a, int b) { return m2[static_cast<size_t>(a) * n + b]; }
    Cost M1(int a, int b) const { return m1[static_cast<size_t>(a) * n + b]; }
    Cost M2(int a, int b) const { return m2[static_cast<size_t>(a) * n + b]; }
};

Graph random_graph(std::mt19937& rng) {
    Graph g;
    g.layers = 3 + static_cast<int>(rng() % 6);  // final position 3..8
    g.at.resize(g.layers + 1);
    for (int p = 0; p <= g.layers; ++p) {
        int width = 1 + static_cast<int>(rng() % 3);
        for (int k = 0; k < width; ++k) {
            int node = static_cast<int>(g.pos.size());
            g.pos.push_back(p);
            g.at[p].push_back(node);
        }
    }
    g.out.resize(g.nodes());
    g.in.resize(g.nodes());
    std::uniform_int_distribution<int> edge_cost(-7, 9);
    for (int p = 0; p < g.layers; ++p)
        for (int from : g.at[p])
            for (int to : g.at[p + 1])
                if (rng() % 100 < 62) {
                    Arc arc{from, to, edge_cost(rng)};
                    g.out[from].push_back(arc);
                    g.in[to].push_back(arc);
                }

    const int n = g.nodes();
    std::vector<unsigned char> reachable(static_cast<size_t>(n) * n, 0);
    for (int a = 0; a < n; ++a) {
        reachable[static_cast<size_t>(a) * n + a] = 1;
        for (int p = g.pos[a]; p < g.layers; ++p)
            for (int from : g.at[p]) {
                if (!reachable[static_cast<size_t>(a) * n + from]) continue;
                for (const Arc& arc : g.out[from])
                    reachable[static_cast<size_t>(a) * n + arc.to] = 1;
            }
    }
    g.branch.assign(static_cast<size_t>(n) * n, INF);
    std::uniform_int_distribution<int> branch_cost(-25, 35);
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b)
            if (g.pos[a] < g.pos[b] && reachable[static_cast<size_t>(a) * n + b] &&
                rng() % 100 < 68)
                g.branch[static_cast<size_t>(a) * n + b] = branch_cost(rng);
    return g;
}

std::vector<Cost> unpaired_paths(const Graph& g) {
    const int n = g.nodes();
    std::vector<Cost> u(static_cast<size_t>(n) * n, INF);
    for (int a = 0; a < n; ++a) u[static_cast<size_t>(a) * n + a] = 0;
    for (int gap = 1; gap <= g.layers; ++gap)
        for (int p = 0; p + gap <= g.layers; ++p)
            for (int a : g.at[p])
                for (int b : g.at[p + gap])
                    for (const Arc& arc : g.out[a])
                        u[static_cast<size_t>(a) * n + b] = std::min(
                            u[static_cast<size_t>(a) * n + b],
                            add(arc.cost, u[static_cast<size_t>(arc.to) * n + b]));
    return u;
}

Tables dense_right(const Graph& g) {
    Tables t(g.nodes());
    for (int gap = 1; gap <= g.layers; ++gap)
        for (int p = 0; p + gap <= g.layers; ++p)
            for (int a : g.at[p])
                for (int b : g.at[p + gap]) {
                    Cost m2 = INF;
                    for (const Arc& arc : g.in[b])
                        m2 = std::min(m2, add(t.M2(a, arc.from), arc.cost));
                    for (int q = p + 1; q < p + gap; ++q)
                        for (int c : g.at[q])
                            m2 = std::min(m2, add(t.M1(a, c), g.B(c, b)));
                    t.M2(a, b) = m2;

                    Cost best = std::min(m2, g.B(a, b));
                    for (const Arc& arc : g.out[a])
                        best = std::min(best, add(arc.cost, t.M1(arc.to, b)));
                    for (const Arc& arc : g.in[b])
                        best = std::min(best, add(t.M1(a, arc.from), arc.cost));
                    t.M1(a, b) = best;
                }
    return t;
}

Tables dense_left(const Graph& g, const std::vector<Cost>& unpaired) {
    const int n = g.nodes();
    Tables t(n);
    for (int gap = 1; gap <= g.layers; ++gap)
        for (int p = 0; p + gap <= g.layers; ++p)
            for (int a : g.at[p])
                for (int b : g.at[p + gap]) {
                    Cost m1 = INF, m2 = INF;
                    for (const Arc& arc : g.out[a]) {
                        m1 = std::min(m1, add(arc.cost, t.M1(arc.to, b)));
                        m2 = std::min(m2, add(arc.cost, t.M2(arc.to, b)));
                    }
                    for (int q = p + 1; q <= p + gap; ++q)
                        for (int c : g.at[q]) {
                            Cost direct = g.B(a, c);
                            m1 = std::min(m1, add(direct, unpaired[static_cast<size_t>(c) * n + b]));
                            m1 = std::min(m1, add(direct, t.M1(c, b)));
                            m2 = std::min(m2, add(direct, t.M1(c, b)));
                        }
                    t.M1(a, b) = m1;
                    t.M2(a, b) = m2;
                }
    return t;
}

struct Candidate { int start; Cost cost; };

Tables sparse_right(const Graph& g, std::int64_t& candidates, std::int64_t& discarded_ties) {
    Tables t(g.nodes());
    std::vector<std::vector<Candidate>> by_end(g.nodes());
    for (int gap = 1; gap <= g.layers; ++gap)
        for (int p = 0; p + gap <= g.layers; ++p)
            for (int a : g.at[p])
                for (int b : g.at[p + gap]) {
                    Cost m2 = INF;
                    for (const Arc& arc : g.in[b])
                        m2 = std::min(m2, add(t.M2(a, arc.from), arc.cost));
                    for (const Candidate& candidate : by_end[b])
                        if (g.pos[candidate.start] > g.pos[a])
                            m2 = std::min(m2, add(t.M1(a, candidate.start), candidate.cost));
                    t.M2(a, b) = m2;

                    Cost alternative = m2;
                    for (const Arc& arc : g.out[a])
                        alternative = std::min(alternative, add(arc.cost, t.M1(arc.to, b)));
                    for (const Arc& arc : g.in[b])
                        alternative = std::min(alternative, add(t.M1(a, arc.from), arc.cost));
                    Cost direct = g.B(a, b);
                    t.M1(a, b) = std::min(alternative, direct);
                    if (direct < alternative) {
                        by_end[b].push_back(Candidate{a, direct});
                        ++candidates;
                    } else if (direct < INF && direct == alternative) {
                        ++discarded_ties;
                    }
                }
    return t;
}

// Independent augmented-DAG shortest path: a transition is either one ordinary unpaired arc or
// one arbitrary direct-branch hyperedge. Branch count is capped at two.
std::pair<Cost, Cost> exhaustive(const Graph& g, int start, int target) {
    const int n = g.nodes();
    std::vector<Cost> distance(static_cast<size_t>(n) * 3, INF);
    distance[static_cast<size_t>(start) * 3] = 0;
    for (int p = g.pos[start]; p < g.pos[target]; ++p)
        for (int node : g.at[p])
            for (int count = 0; count <= 2; ++count) {
                Cost here = distance[static_cast<size_t>(node) * 3 + count];
                if (here >= INF) continue;
                for (const Arc& arc : g.out[node])
                    if (g.pos[arc.to] <= g.pos[target])
                        distance[static_cast<size_t>(arc.to) * 3 + count] = std::min(
                            distance[static_cast<size_t>(arc.to) * 3 + count], add(here, arc.cost));
                for (int q = p + 1; q <= g.pos[target]; ++q)
                    for (int end : g.at[q]) {
                        Cost branch = g.B(node, end);
                        if (branch >= INF) continue;
                        int next_count = std::min(2, count + 1);
                        distance[static_cast<size_t>(end) * 3 + next_count] = std::min(
                            distance[static_cast<size_t>(end) * 3 + next_count], add(here, branch));
                    }
            }
    Cost one = std::min(distance[static_cast<size_t>(target) * 3 + 1],
                        distance[static_cast<size_t>(target) * 3 + 2]);
    return {one, distance[static_cast<size_t>(target) * 3 + 2]};
}

}  // namespace

int main() {
    std::mt19937 rng(0x5a17eU);
    int failures = 0;
    std::int64_t cells = 0, candidates = 0, ties = 0;
    constexpr int kGraphs = 500;
    for (int trial = 0; trial < kGraphs; ++trial) {
        Graph graph = random_graph(rng);
        std::vector<Cost> unpaired = unpaired_paths(graph);
        Tables left = dense_left(graph, unpaired);
        Tables right = dense_right(graph);
        Tables sparse = sparse_right(graph, candidates, ties);
        for (int p = 0; p < graph.layers; ++p)
            for (int q = p + 1; q <= graph.layers; ++q)
                for (int a : graph.at[p])
                    for (int b : graph.at[q]) {
                        auto exact = exhaustive(graph, a, b);
                        bool ok = left.M1(a, b) == exact.first && left.M2(a, b) == exact.second &&
                                  right.M1(a, b) == exact.first && right.M2(a, b) == exact.second &&
                                  sparse.M1(a, b) == exact.first && sparse.M2(a, b) == exact.second;
                        if (!ok) {
                            ++failures;
                            if (failures <= 5)
                                std::printf("FAIL trial=%d cell=(%d,%d) exact=(%lld,%lld) "
                                            "left=(%lld,%lld) right=(%lld,%lld) sparse=(%lld,%lld)\n",
                                            trial, a, b, (long long)exact.first,
                                            (long long)exact.second, (long long)left.M1(a, b),
                                            (long long)left.M2(a, b), (long long)right.M1(a, b),
                                            (long long)right.M2(a, b), (long long)sparse.M1(a, b),
                                            (long long)sparse.M2(a, b));
                        }
                        ++cells;
                    }
    }
    std::printf("abstract sparse kernel: graphs=%d cells=%lld candidates=%lld discarded_ties=%lld "
                "failures=%d\n", kGraphs, (long long)cells, (long long)candidates,
                (long long)ties, failures);
    return failures == 0 ? 0 : 1;
}
