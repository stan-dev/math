//=======================================================================
// Copyright (c) 2024 Andrea Cassioli
// Author: Andrea Cassioli <cassioliandre@gmail.com>
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/r_c_shortest_paths.hpp>

#include <boost/core/lightweight_test.hpp>

struct edge_prop {
  edge_prop(int c, int t) : cost(c), time(t) {}

  int cost;
  int time;
};

using graph_type = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, std::string, edge_prop>;
using vertex_type = graph_type::vertex_descriptor;
using edge_type = graph_type::edge_descriptor;
using path_type = std::vector<edge_type>;

class dominance
{
 public:
    template<class C>
    inline bool operator()(const C& res_cont_1, const C& res_cont_2) const {
      return res_cont_1.cost <= res_cont_2.cost && res_cont_1.time <= res_cont_2.time;
    }
};

struct resource_container
{
  resource_container(int c, int t) : cost(c), time(t) {}

  int cost; // minimise cost
  int time; // a fake time constraints
};


bool operator==(
        const resource_container& res_cont_1, const resource_container& res_cont_2)
{
    return (res_cont_1.cost == res_cont_2.cost
            && res_cont_1.time == res_cont_2.time);
}

bool operator<(
        const resource_container& res_cont_1, const resource_container& res_cont_2)
{
    if (res_cont_1.cost > res_cont_2.cost)
        return false;
    if (res_cont_1.cost == res_cont_2.cost)
        return res_cont_1.time < res_cont_2.time;
    return true;
}

struct extension_function {
  template <typename GraphType>
  bool operator()(const GraphType& graph,
                  resource_container& new_cont,
                  const resource_container& old_cont,
                  const typename GraphType::edge_descriptor& edge) const {
    new_cont = old_cont;
    new_cont.cost = old_cont.cost + graph[edge].cost;
    // here I could check tme constraint, but for this example does not matter
    new_cont.time = old_cont.time + graph[edge].time;
    return true;
  }
};

resource_container run_rcsp(const graph_type& graph, vertex_type source, vertex_type target) {
  const auto vertex_index_map = boost::get(boost::vertex_index, graph);
  const auto edge_index_map = boost::get(boost::edge_all, graph);
  boost::default_r_c_shortest_paths_allocator label_allocator{};

  path_type single_solution;
  resource_container single_resource(0, 0);
  const resource_container start_resource(0, 0);
  boost::r_c_shortest_paths(graph, vertex_index_map, edge_index_map, source, target, single_solution, single_resource,
                            start_resource, extension_function{}, dominance{}, label_allocator,
                            boost::default_r_c_shortest_paths_visitor());

  return single_resource;
}

int main() {
  graph_type graph;

  /*
        (1,0)       (10, 1)
       /-----> [A] ------\
    [s]                   [t]
       \-----> [B] ------/
        (2, 1)      (3,1)

    The shortest path is s->B->t with cost 5.
  */

  const auto s = boost::add_vertex("s", graph);
  const auto A = boost::add_vertex("A", graph);
  const auto B = boost::add_vertex("B", graph);
  const auto t = boost::add_vertex("t", graph);

  boost::add_edge(s, A, edge_prop(1, 0), graph);
  boost::add_edge(s, B, edge_prop(2, 1), graph);
  boost::add_edge(A, t, edge_prop(10, 0), graph);
  boost::add_edge(B, t, edge_prop(3, 1), graph);

  const auto solution_resource_container = run_rcsp(graph, s, t);

  BOOST_TEST(solution_resource_container.cost == 5);
  return boost::report_errors();
}
