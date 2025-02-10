// Copyright 2004-5 Trustees of Indiana University
//

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//
// graphviz_test.cpp - Test cases for the Graphviz DOT Language reader
//

// Author: Ronald Garcia

#define BOOST_GRAPHVIZ_USE_ISTREAM
#include <boost/graph/graphviz.hpp>
#include <boost/assign/std/map.hpp>
#include <boost/assign.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <iostream>
#include <cmath>

typedef std::string node_t;
typedef std::pair< node_t, node_t > edge_t;

typedef float Mass;
typedef double Weight;
typedef std::map< node_t, Mass > expected_masses_t;
typedef std::map< edge_t, Weight > expected_weights_t;
#define MAP_MASSES boost::assign::list_of< std::pair< node_t, Mass > >
#define MAP_WEIGHTS boost::assign::list_of< std::pair< edge_t, Weight > >

struct Fixture
{
    std::string graphviz_text;
    size_t correct_num_vertices;
    expected_masses_t masses;
    expected_weights_t weights;
};

namespace Samples
{
namespace Directed
{
    static Fixture const basic {
        "digraph { a  node [mass = 7.7] c e [mass = 6.66] }",
        3,
        MAP_MASSES("a", 0.0f)("c", 7.7f)("e", 6.66f),
        expected_weights_t(),
    };

    static Fixture const basic_aliased {
        "digraph { a  node [mass = 7.7] \"a\" e [mass = 6.66] }",
        2,
        MAP_MASSES("a", 0.0f)("e", 6.66f),
        expected_weights_t(),
    };

    static Fixture const full {
        "digraph { a -> b eDge [weight = 7.7] "
        "c -> d e-> f [weight = 6.66] "
        "d ->e->a [weight=.5]}",
        6,
        expected_masses_t(),
        MAP_WEIGHTS(edge_t("a", "b"), 0.0)(edge_t("c", "d"), 7.7)(
            edge_t("e", "f"), 6.66)(edge_t("d", "e"), 0.5)(
            edge_t("e", "a"), 0.5),
    };
}

namespace Undirected
{
    static Fixture const basic {
        "graph { a  nodE [mass = 7.7] c e [mass =\\\n6.66] }",
        3,
        MAP_MASSES("a", 0.0f)("c", 7.7f)("e", 6.66f),
        expected_weights_t(),
    };

    static Fixture const full {
        "graph { a -- b eDge [weight = 7.7] "
        "c -- d e -- f [weight = 6.66] }",
        6,
        expected_masses_t(),
        MAP_WEIGHTS(edge_t("a", "b"), 0.0)(edge_t("c", "d"), 7.7)(
            edge_t("e", "f"), 6.66),
    };
}

Fixture const all_directed[] {
    Directed::basic,
    Directed::basic_aliased,
    Directed::full,
};
Fixture const all_undirected[] {
    Undirected::basic,
    Undirected::full,
};
} // Samples

namespace ComparisonDriver
{
template < class T > class close_to
{
public:
    explicit close_to(T f) : f_(f) { assert(f >= 0); }

    bool operator()(T l, T r) const
    {
        using std::abs;
        return abs(l - r) <= f_ * (std::max)(abs(l), abs(r));
    }

private:
    T f_;
};

template < typename graph_t, typename NameMap, typename MassMap,
    typename WeightMap >
bool test_graph(std::string const& text, graph_t& graph,
    std::size_t correct_num_vertices, expected_masses_t const& expected_masses,
    expected_weights_t const& expected_weights, std::string const& node_id,
    std::string const& g_name, NameMap name_map, MassMap mass_map,
    WeightMap weight_map)
{
    // Construct a graph and set up the dynamic_property_maps.
    boost::dynamic_properties dp(boost::ignore_other_properties);
    dp.property(node_id, name_map);
    dp.property("mass", mass_map);
    dp.property("weight", weight_map);

    boost::ref_property_map< graph_t*, std::string > gname(
        get_property(graph, boost::graph_name));
    dp.property("name", gname);

    bool result = true;
#ifdef BOOST_GRAPHVIZ_USE_ISTREAM
    std::istringstream is(text);
    if (read_graphviz(is, graph, dp, node_id))
#else
    if (read_graphviz(text.begin(), text.end(), graph, dp, node_id))
#endif
    {
        // check correct vertex count
        BOOST_TEST_EQ(num_vertices(graph), correct_num_vertices);
        // check masses
        if (!expected_masses.empty())
        {
            // assume that all the masses have been set
            // for each vertex:
            typename boost::graph_traits< graph_t >::vertex_iterator i, j;
            for (boost::tie(i, j) = vertices(graph); i != j; ++i)
            {
                std::string node_name = get(name_map, *i);
                Mass node_mass = get(mass_map, *i);
                BOOST_TEST(
                    expected_masses.find(node_name) != expected_masses.end());
                Mass ref_mass = expected_masses.find(node_name)->second;
                //  - compare the mass to the result in the table
                BOOST_TEST_WITH(node_mass, ref_mass, close_to< Mass >(0.0001f));
            }
        }
        // check weights
        if (!expected_weights.empty())
        {
            // assume that all weights have been set
            /// for each edge:
            typename boost::graph_traits< graph_t >::edge_iterator i, j;
            for (boost::tie(i, j) = edges(graph); i != j; ++i)
            {
                //  - get its name
                std::pair< std::string, std::string > edge_name
                    = make_pair(get(name_map, source(*i, graph)),
                        get(name_map, target(*i, graph)));
                // - get its weight
                Weight edge_weight = get(weight_map, *i);
                BOOST_TEST(
                    expected_weights.find(edge_name) != expected_weights.end());
                Weight ref_weight = expected_weights.find(edge_name)->second;
                // - compare the weight to the result in the table
                BOOST_TEST_WITH(
                    edge_weight, ref_weight, close_to< Weight >(0.000001));
            }
        }
        if (!g_name.empty())
        {
            std::string parsed_name = get_property(graph, boost::graph_name);
            BOOST_TEST(parsed_name == g_name);
        }
    }
    else
    {
        std::cerr << "Parsing Failed!\n";
        result = false;
    }

    return result;
}

template < typename graph_t, typename NameMap, typename MassMap,
    typename WeightMap >
bool test_graph(Fixture const& sample, graph_t& g, std::string const& node_id,
    std::string const& g_name, NameMap name_map, MassMap mass_map,
    WeightMap weight_map)
{
    return test_graph(sample.graphviz_text, g, sample.correct_num_vertices,
        sample.masses, sample.weights, node_id, g_name, name_map, mass_map,
        weight_map);
}

template < typename graph_t >
bool test_graph(std::string const& dottext, std::size_t correct_num_vertices,
    expected_masses_t const& masses, expected_weights_t const& weights,
    std::string const& node_id = "node_id", std::string const& g_name = "")
{
    graph_t g;
    return test_graph(dottext, g, correct_num_vertices, masses, weights,
        node_id, g_name, get(boost::vertex_name, g),
        get(boost::vertex_color, g), get(boost::edge_weight, g));
}

template < typename graph_t >
bool test_graph(Fixture const& sample, std::string const& node_id = "node_id",
    std::string const& g_name = "")
{
    graph_t g;
    return test_graph(sample.graphviz_text, g, sample.correct_num_vertices,
        sample.masses, sample.weights, node_id, g_name, //
        get(boost::vertex_name, g), //
        get(boost::vertex_color, g), //
        get(boost::edge_weight, g) //
    );
}
}

namespace Models
{
typedef boost::property< boost::vertex_name_t, std::string,
    boost::property< boost::vertex_color_t, Mass > >
    vertex_p;
typedef boost::property< boost::edge_weight_t, Weight > edge_p;
typedef boost::property< boost::graph_name_t, std::string > graph_p;

typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS,
    vertex_p, edge_p, graph_p >
    Graph;
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS,
    vertex_p, edge_p, graph_p >
    DiGraph;

typedef boost::adjacency_list< boost::setS, boost::vecS, boost::directedS,
    vertex_p, edge_p, graph_p >
    DiGraphNoParallel;

struct VertexBundle
{
    std::string name;
    Mass mass;
};
struct EdgeBundle
{
    Weight weight;
};

typedef boost::compressed_sparse_row_graph< boost::directedS, VertexBundle,
    EdgeBundle, graph_p >
    CSRBundledGraph;
typedef boost::compressed_sparse_row_graph< boost::directedS,
    boost::no_property, boost::no_property, graph_p >
    CSRGraph;
} // Models

// SEHE I intended for this to be able to pass __FUNCTION__ for reporting
// failures, but there seems to be no way to achieve this with lightweight_test
#define TEST_GRAPH(Model, ...) \
    BOOST_TEST((ComparisonDriver::test_graph< Model >(__VA_ARGS__)));

// Basic directed graph tests
void test_basic_directed_graph_1()
{
    TEST_GRAPH(Models::DiGraph, Samples::Directed::basic);
}

void test_basic_directed_graph_2()
{
    TEST_GRAPH(Models::DiGraph, Samples::Directed::basic_aliased);
}

void test_basic_directed_graph_3()
{
    TEST_GRAPH(Models::DiGraph, Samples::Directed::full);
}

// undirected graph with alternate node_id property name
void test_undirected_graph_alternate_node_id()
{
    TEST_GRAPH(Models::Graph, Samples::Undirected::basic, "nodenames");
}

// Basic undirected graph tests
void test_basic_undirected_graph_1()
{
    TEST_GRAPH(Models::Graph, Samples::Undirected::basic);
}

void test_basic_undirected_graph_2()
{
    TEST_GRAPH(Models::Graph, Samples::Undirected::full);
}

// Mismatch directed graph test
void test_mismatch_directed_graph()
{
    BOOST_TEST_THROWS(TEST_GRAPH(Models::DiGraph, Samples::Undirected::basic),
        boost::undirected_graph_error);
}

// Mismatch undirected graph test
void test_mismatch_undirected_graph()
{
    BOOST_TEST_THROWS(TEST_GRAPH(Models::Graph, Samples::Directed::basic),
        boost::directed_graph_error);
}

// Complain about parallel edges
void test_parallel_edges()
{
    Fixture parallel {
        "diGraph { a -> b [weight = 7.7]  a -> b [weight = 7.7] }",
        2,
        expected_masses_t(),
        MAP_WEIGHTS(edge_t("a", "b"), 7.7),
    };
    TEST_GRAPH(Models::DiGraph, parallel);
    BOOST_TEST_THROWS(TEST_GRAPH(Models::DiGraphNoParallel, parallel),
        boost::bad_parallel_edge);
}

// Graph Property Test 1
void test_graph_property_test_1()
{
    Fixture named {
        "digraph { graph [name=\"foo \\\"escaped\\\"\"]  a  c e [mass = 6.66] "
        "}",
        3,
        MAP_MASSES("a", 0.0f)("c", 0.0f)("e", 6.66f),
        expected_weights_t(),
    };
    TEST_GRAPH(Models::DiGraph, named, "", "foo \"escaped\"");
}

// Graph Property Test 2
void test_graph_property_test_2()
{
    Fixture named {
        "digraph { name=\"fo\"+ \"\\\no\"  a  c e [mass = 6.66] }",
        3,
        MAP_MASSES("a", 0.0f)("c", 0.0f)("e", 6.66f),
        expected_weights_t(),
    };
    TEST_GRAPH(Models::DiGraph, named, "", "foo"); // SEHE why not "foo\no"?
}

// Graph Property Test 3 (HTML)
void test_graph_property_test_3()
{
    std::string graph_name
        = "<html title=\"x'\" title2='y\"'>foo<b><![CDATA[><bad "
          "tag&>]]>bar</b>\n<br/>\nbaz</html>";
    Fixture html_named {
        "digraph { name=" + graph_name + "  a  c e [mass = 6.66] }",
        3,
        MAP_MASSES("a", 0.0f)("c", 0.0f)("e", 6.66f),
        expected_weights_t(),
    };
    TEST_GRAPH(Models::DiGraph, html_named, "", graph_name);
}

// Comments embedded in strings
void test_comments_embedded_in_strings()
{
    std::string gv("digraph { "
                   "a0 [ label = \"//depot/path/to/file_14#4\" ];"
                   "a1 [ label = \"//depot/path/to/file_29#9\" ];"
                   "a0 -> a1 [ color=gray ];"
                   "}");
    TEST_GRAPH(
        Models::DiGraph, gv, 2, expected_masses_t(), expected_weights_t());
}

void test_basic_csr_directed_graph()
{
    auto sample = Samples::Directed::full;

    typedef Models::CSRBundledGraph graph_t;
    graph_t g;
#ifdef UNSTABLE_PROPERTY_MAPS_FIXED // https://github.com/boostorg/graph/issues/373
    BOOST_TEST((test_graph(gs, g, 6, mass_map_t(), weights, "node_id", "",
        get(&vertex_p_bundled::name, g), // Warning, currently broken
        get(&vertex_p_bundled::color, g), // Warning, currently broken
        get(&edge_p_bundled::weight, g)) // Warning, currently broken
        ));
#else
    typedef graph_t::vertex_descriptor V;
    typedef graph_t::edge_descriptor E;
    TEST_GRAPH(graph_t, sample, g, "node_id", "", //
        boost::make_function_property_map< V >(
            [&g](V v) -> std::string& { return g[v].name; }),
        boost::make_function_property_map< V >(
            [&g](V v) -> Mass& { return g[v].mass; }),
        boost::make_function_property_map< E >(
            [&g](E e) -> Weight& { return g[e].weight; }) //
    );
#endif
}

void test_basic_csr_directed_graph_ext_props()
{
    auto sample = Samples::Directed::full;
    using Models::CSRGraph;
    CSRGraph g;
    boost::property_map< CSRGraph, boost::vertex_index_t >::const_type vidx
        = get(boost::vertex_index, g);
    boost::property_map< CSRGraph, boost::edge_index_t >::const_type eidx
        = get(boost::edge_index, g);

    boost::vector_property_map< std::string,
        boost::property_map< CSRGraph, boost::vertex_index_t >::const_type >
        vertex_name(vidx);
    boost::vector_property_map< Mass,
        boost::property_map< CSRGraph, boost::vertex_index_t >::const_type >
        vertex_mass(vidx);
    boost::vector_property_map< Weight,
        boost::property_map< CSRGraph, boost::edge_index_t >::const_type >
        edge_weight(eidx);

    TEST_GRAPH(CSRGraph, sample, g, "node_id", "", vertex_name, vertex_mass,
        edge_weight);
}

void test_subgraphs()
{
    // on the BGL side, the new parser doesn't support subgraphs
    // however, the docs promise to support reading them on the input side as
    // "syntactic sugar".
    for (auto gv : {
             Fixture { "digraph {}" },
             Fixture { "digraph { 1 -> {} }", 1 },
             Fixture { "digraph { 1 -> {2} }", 2 },
             Fixture { "digraph { 1; { 2; 3; } }", 3 },
             Fixture { "digraph { { 2; 3; } 1; }", 3 },
             Fixture { "digraph { 1; subgraph { 2; 3; } }", 3 },
             Fixture { "digraph { 1 -> subgraph { 2; 3; } }", 3 },
             Fixture { "digraph { 1 -> subgraph hello { 2; 3; } }", 3 },
             Fixture { "digraph { 1 -> subgraph clust_Hello { 2; 3; } }", 3 },
             Fixture { "digraph { 1 -> subgraph \"hello\" { 2; 3; } }", 3 },
             Fixture {
                 "digraph { {2} -> subgraph \"hello\" {{{{ 2; 3; }}}} }", 2 },
         })
    {
        TEST_GRAPH(Models::DiGraph, gv);
    }
}

void test_subgraph_nesting_limit() // issue #364
{
    auto independent_nests = [=](unsigned level)
    {
        auto sg = std::string(level, '{') + " 2; 3; " + std::string(level, '}');
        ComparisonDriver::test_graph< Models::DiGraph >(
            { "digraph{1->" + sg + "}", 3 });
        ComparisonDriver::test_graph< Models::DiGraph >(
            { "digraph{1->" + sg + ";4->" + sg + "}", 4 });
    };

    constexpr unsigned limit = 255;
    independent_nests(1);
    independent_nests(limit / 2);
    independent_nests(limit - 1);
    independent_nests(limit); // edge-case
    BOOST_TEST_THROWS(independent_nests(limit + 1), boost::bad_graphviz_syntax);
}

int main()
{
    test_basic_directed_graph_1();
    test_basic_directed_graph_2();
    test_basic_directed_graph_3();
    test_undirected_graph_alternate_node_id();
    test_basic_undirected_graph_1();
    test_basic_undirected_graph_2();
    test_mismatch_directed_graph();
    test_mismatch_undirected_graph();
    test_parallel_edges();
    test_graph_property_test_1();
    test_graph_property_test_2();
    test_graph_property_test_3();
    test_comments_embedded_in_strings();
    test_basic_csr_directed_graph_ext_props();
    test_basic_csr_directed_graph();
    test_subgraphs();
    test_subgraph_nesting_limit();
    return boost::report_errors();
}
