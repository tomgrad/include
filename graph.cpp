#include "graph.h"

#include <algorithm>
#include <fstream>

void Graph::add_edge(const size_t i, const size_t j)
{
    A[i].push_back(j);
    A[j].push_back(i);
}

void Graph::add_directed_edge(const size_t from, const size_t to)
{
    A[from].push_back(to);
}

void Graph::reset(const size_t size = 0)
{
    A.clear();
    A.resize(size);
}

bool Graph::is_connected(const size_t i, const size_t j) const
{
    return std::find(begin(A[i]), end(A[i]), j) != end(A[i]);
}

std::ostream &operator<<(std::ostream &os, const Graph &G)
{
    unsigned n = 0;
    for (auto &a : G.A)
    {
        os << n++ << ": ";
        for (auto i : a)
            os << i << " ";
        os << '\n';
    }

    return os;
}

namespace graph
{
std::vector<unsigned> get_k(const Graph &G)
{
    std::vector<unsigned> result;
    result.reserve(G.A.size());

    for (const auto &a : G.A)
        result.push_back(a.size());
    return result;
}

size_t get_k(const Graph &G, const size_t i) { return G.A[i].size(); }
} // namespace graph

std::vector<unsigned> gen_discrete_powerlaw(unsigned xmin, double gamma,
                                            size_t N, URNG &rng)
{

    // auto Pc = [](auto x, auto xmin, auto g) { return pow(x / xmin, -g + 1); };

    // auto P = [&Pc](auto x, auto xmin, auto g) {
    //   return (g - 1) * pow(xmin, g - 1) * pow(x, -g);
    //   // return Pc(x, xmin, g) - Pc(x + 1, xmin, g);
    // };

    auto P = [](auto x, auto g) { return pow(x, -g); };
    std::vector<double> Pv(N);
    for (size_t i = 0; i < N; ++i)
        Pv[i] = P(i, gamma);
    auto A = std::accumulate(std::begin(Pv) + xmin, std::end(Pv), 0.0);
    for (auto &x : Pv)
        x /= A;

    auto prob = std::uniform_real_distribution<>(0, 1);
    std::vector<unsigned> v(N);

    for (size_t i = 0; i < N; ++i)
    {
        unsigned x = xmin;
        double sum = Pv[x];
        auto r = prob(rng);
        while (sum < r)
            sum += Pv[++x];
        v[i] = x;
    }
    return v;
}

namespace graph
{

/*  Erdos-Renyi  */
void make_er(Graph &G, const double p, URNG &rng)
{
    auto N = G.A.size();
    auto prob = std::uniform_real_distribution<>(0, 1);

    for (size_t i = 0; i != N - 1; ++i)
        for (size_t j = i + 1; j != N; ++j)
            if (prob(rng) < p)
                G.add_edge(i, j);
}

/*  Barabasi-Albert  */
void make_ba(Graph &G, const unsigned m, const unsigned m0, URNG &rng)
{
    if (m0 < 1 || m > m0)
        return;
    auto N = G.A.size();
    auto dist = std::uniform_int_distribution<size_t>();
    auto randint = [&](auto range) { return dist(rng) % range; };

    auto is_in = [](auto val, auto cont) {
        return find(begin(cont), end(cont), val) != end(cont);
    };

    for (auto i = 0u; i < m0 - 1; ++i)
        for (auto j = i + 1; j < m0; ++j)
            G.add_edge(i, j);
    unsigned sumk = m0 * (m0 - 1);
    for (size_t i = m0; i < N; ++i)
    {
        std::vector<size_t> nodes;
        while (nodes.size() != m)
        {
            int r = randint(sumk); // unsigned

            size_t j;
            while (r > 0)
            {
                j = randint(i);
                r -= graph::get_k(G, j);
            }
            if (!is_in(j, nodes))
                nodes.push_back(j);
        }
        for (auto j : nodes)
            G.add_edge(i, j);
        sumk += m << 1;
    }
}

/*  configuration model */
void make_conf(Graph &G, std::vector<unsigned> &deg, URNG &rng)
{

    if (accumulate(begin(deg), end(deg), 0) % 2)
        ++deg[0];

    std::vector<unsigned> stub_list;
    for (auto i = 0u; i != deg.size(); ++i)
        for (auto k = 0u; k != deg[i]; ++k)
            stub_list.push_back(i);

    auto rand_int = std::uniform_int_distribution<unsigned>();
    auto draw_int = [&](size_t x) { return rand_int(rng) % x; };

    std::random_shuffle(stub_list.begin(), stub_list.end(), draw_int);

    if (stub_list.size() % 2)
    {
        stub_list.push_back(draw_int(deg.size()));
    }

    G.reset(deg.size());
    while (!stub_list.empty())
    {
        auto i = *(stub_list.end() - 1);
        auto j = *(stub_list.end() - 2);

        G.add_edge(i, j);
        stub_list.pop_back();
        stub_list.pop_back();
    }
}

// TODO: stopniowe i progresywne usuwanie zlych linkow
void fix(Graph &G, unsigned iter, URNG &rng)
{
    unsigned n_loops = 0, n_dupes = 0;
    vert_vec free_nodes;
    auto rand_int = std::uniform_int_distribution<unsigned>();
    auto draw_int = [&](size_t x) { return rand_int(rng) % x; };

    auto remove_it = [](auto &v, auto &it) {
        std::swap(*it, v.back());
        v.pop_back();
    };

    for (unsigned i = 0; i < G.A.size(); ++i)
    {
        auto &a = G.A[i];

        // randomize loops
        auto it = begin(a);
        do
        {
            it = std::find(begin(a), end(a), i);

            if (it != a.end()) // jest petla
            {
                ++n_loops;
                remove_it(a, it);
                it = std::find(begin(a), end(a), i);
                remove_it(a, it);

                if (iter > 1)
                {
                    unsigned j = 0, k = 0;
                    do
                    {
                        j = draw_int(G.A.size());
                        k = G.A[j][draw_int(G.A[j].size())];
                    } while (i == j || i == k || j == k);
                    auto &aj = G.A[j];
                    auto &ak = G.A[k];

                    auto jt = std::find(begin(aj), end(aj), k);
                    auto kt = std::find(begin(ak), end(ak), j);

                    remove_it(aj, jt);
                    remove_it(ak, kt);

                    G.add_edge(i, j);
                    G.add_edge(i, k);
                    //   std::cout << i << " " << j << " " << k << "\n";
                }
            }

        } while (it != a.end());
    }

    for (unsigned i = 0; i < G.A.size(); ++i)
    {
        auto &a = G.A[i];

        // find multiedges
        //   auto it = find_first_of(begin(a), end(a),begin(a), end(a));
        std::sort(begin(a), end(a));

        auto it = begin(a);

        do
        {

            it = std::adjacent_find(begin(a), end(a));
            if (it != a.end()) // jest duplikat
            {
                ++n_dupes;
                auto j = *it;
                auto &b = G.A[j];
                auto jt = std::find(begin(b), end(b), i);

                remove_it(a, it);
                remove_it(b, jt);
                if (iter > 1)
                {
                    free_nodes.push_back(i);
                    free_nodes.push_back(j);
                }
            }
        } while (it != end(a));
    }

    // std::cout << "loops: " << n_loops << " dupes: " << n_dupes << '\n';
    if (!n_loops && !n_dupes)
        return;

    // remove some good links

    for (unsigned i = 0; i < free_nodes.size() / 10 + 1; ++i)
    {
        auto j = draw_int(G.A.size());
        if (!G.A[j].empty())
        {
            auto k = G.A[j][draw_int(G.A[j].size())];
            auto &aj = G.A[j];
            auto &ak = G.A[k];
            auto jt = std::find(begin(aj), end(aj), k);
            auto kt = std::find(begin(ak), end(ak), j);

            remove_it(aj, jt);
            remove_it(ak, kt);
            free_nodes.push_back(j);
            free_nodes.push_back(k);
        }
    }

    // std::cout << "Free nodes: ";
    // for (auto i : free_nodes)
    //   std::cout << i << " ";

    std::random_shuffle(free_nodes.begin(), free_nodes.end(), draw_int);

    for (unsigned i = 0; i < free_nodes.size(); i += 2)
        G.add_edge(free_nodes[i], free_nodes[i + 1]);

    if (iter)
    {
        graph::fix(G, --iter, rng);
    }
}

/*  Random regular graph  */
void make_rr(Graph &G, const unsigned k, URNG &rng)
{

    auto N = G.A.size();
    std::vector<unsigned> v(N, k);

    graph::make_conf(G, v, rng);

    graph::fix(G, N, rng);
}

/*  Scale free graph (conf. model)  */
void make_sf(Graph &G, const double gamma, const unsigned kmin, URNG &rng)
{
    auto N = G.A.size();
    auto v = gen_discrete_powerlaw(kmin, gamma, N, rng);
    sort(rbegin(v), rend(v));
    graph::make_conf(G, v, rng);
    graph::fix(G, N, rng);
}

} // namespace graph

/*
std::istream &operator>>(std::istream &is, HGraph &G) {

  std::string line;
  std::getline(is, line);
  size_t n;
  std::stringstream(line) >> n;
  hgraph::reset(G, n);
  while (std::getline(is, line)) {
    edge e;
    std::stringstream ss(line);
    while (ss.good()) {
      int i;
      if (ss >> i) {
        e.push_back(i);
      }
    }
    hgraph::add_edge(e, G);
  }
  return is;
}
*/

/*

void makeSqr(Graph &G, int radius = 1, bool periodic = true) {
    int L = sqrt(num_vertices(G));
    int i, j, r;

    if (periodic) {
        for (r = 1; r <= radius; ++r)
            for (i = 0; i < L; i++)
                for (j = 0; j < L; j++) {
                    add_edge(i * L + j, i * L + (j + r) % L, G);
                    add_edge(i * L + j, (i + r) % L * L + j, G);
                    add_edge(i * L + (j + r) % L, i * L + j, G);
                    add_edge((i + r) % L * L + j, i * L + j, G);
                }
    } else {
        for (r = 1; r <= radius; ++r)
            for (i = 0; i < L; i++)
                for (j = 0; j < L; j++) {
                    if (i < L - r) {
                        add_edge(j * L + i, j * L + i + r, G);
                        add_edge(j * L + i + r, j * L + i, G);
                    }
                    if (j < L - r) {
                        add_edge(j * L + i, (j + r) * L + i, G);
                        add_edge((j + r) * L + i, j * L + i, G);
                    }
                }
    }
}
*/
