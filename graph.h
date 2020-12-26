#ifndef _GRAPH_H
#define _GRAPH_H

#include <iostream>
#include <random>
#include <vector>

using URNG = std::mt19937;

typedef std::vector<unsigned> vert_vec;

struct Graph
{
  Graph(const size_t size = 0) : A(size) {}
  void add_edge(const size_t, const size_t);
  void add_directed_edge(const size_t, const size_t);
  void reset(const size_t);
  bool is_connected(const size_t, const size_t) const;
  std::vector<vert_vec> A; // adjacency list
};

std::ostream &operator<<(std::ostream &, const Graph &);

std::vector<unsigned> gen_discrete_powerlaw(unsigned, double, size_t, URNG &);

namespace graph
{
std::vector<unsigned> get_k(const Graph &);
size_t get_k(const Graph &, const size_t);

/*  Erdos-Renyi  */
void make_er(Graph &, const double, URNG &);

/*  Barabasi-Albert  */
void make_ba(Graph &, const unsigned, const unsigned, URNG &);

/*  configuration model */
void make_conf(Graph &, std::vector<unsigned> &, URNG &);

void fix(Graph &, unsigned, URNG &);

/*  Random regular graph  */
void make_rr(Graph &, const unsigned, URNG &);

/*  Scale free graph (conf. model)  */
void make_sf(Graph &, const double, const unsigned, URNG &);

} // namespace graph

#endif // _GRAPH_H
