// Boost.Geometry
// Unit Test

// Copyright (c) 2023 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/geometry.hpp>

#include "helper_functions.hpp"

#include "cnc_container.hpp"
#include "cnc_linestring.hpp"
#include "cnc_multi_linestring.hpp"

#include "adapt_cnc_container.hpp"
#include "adapt_cnc_linestring.hpp"
#include "adapt_cnc_multi_linestring.hpp"

#include <geometry_test_common.hpp>

#include <sstream>

namespace bg = boost::geometry;

template <typename P>
void test_const()
{
    using multi_t = cnc_multi_linestring<cnc_linestring<P>>;
    using point_t = typename bg::point_type<multi_t>::type;

    boost::ignore_unused<point_t>();

    BOOST_CONCEPT_ASSERT( (bg::concepts::ConstMultiLinestring<multi_t>) );

    multi_t geo;
    geo.custom_resize(2);
    fill(geo.custom_get(0), {{0, 0}, {5, 5}, {7, 3}, {9, 5}, {10, 10}});
    fill(geo.custom_get(1), {{0, 10}, {2, 8}, {4, 9}});

#if defined(TEST_FAIL_CNC)
    // This should NOT work.
    auto copy = geo;
#endif

    BOOST_CHECK_EQUAL(8u, bg::num_points(geo));
    BOOST_CHECK_EQUAL(false, bg::is_empty(geo));
    BOOST_CHECK_EQUAL(false, bg::is_convex(geo));
    BOOST_CHECK_EQUAL(true, bg::is_simple(geo));
    BOOST_CHECK_EQUAL(true, bg::is_valid(geo));
    BOOST_CHECK_EQUAL(6u, bg::num_segments(geo));
    BOOST_CHECK_EQUAL(2u, bg::num_geometries(geo));

    // Check floating point properties, they have to compile, matching exactly is done elsewhere
    BOOST_CHECK_GT(bg::length(geo), 10.0);

#if defined(BOOST_GEOMETRY_TODO_SUPPORT_CUSTOM_BUFFER_RESULT)
    using result_t = cnc_multi_polygon<cnc_polygon<P>>;
#else
    using result_t = boost::geometry::model::multi_polygon<boost::geometry::model::polygon<P>>;
#endif

    auto const buffered = buffer<result_t>(geo, 0.5);
    BOOST_CHECK_GT(bg::area(buffered), bg::length(geo));
    BOOST_CHECK_EQUAL(true, bg::within(geo, buffered));

    std::ostringstream out;
    out << "MULTI LS: " << bg::wkt(geo)
        << " centroid: " << bg::wkt(centroid(geo))
        << " envelope: " << bg::wkt(envelope(geo))
        << " hull: " << bg::wkt(hull(geo));

    // Test Boost.Range behavior
    for (auto it = boost::begin(geo); it != boost::end(geo); ++it)
    {
        out << " " << bg::wkt(*it) << " " << bg::dsv(*it);
    }
    out << std::endl;

#if 0
    // Test std::begin/end (does not work)
    for (auto it = std::begin(geo); it != std::end(geo); ++it)
    {
        out << " " << bg::wkt(*it) << " " << bg::dsv(*it);
    }
    out << std::endl;

    // Test range based for loops (does not work)
    for (auto const& point : geo)
    {
        out << " " << bg::wkt(point);
    }
    out << std::endl;
#endif
    auto const test_string = out.str();
    auto const test_length = test_string.length();
    BOOST_CHECK_MESSAGE(test_length > 150,
        "detected: '" << test_string << "' (" << test_length << ")");
}

template <typename P>
void test_mutable()
{
    using multi_t = cnc_multi_linestring<cnc_linestring<P>>;
    using point_t = typename bg::point_type<multi_t>::type;

    boost::ignore_unused<point_t>();

    BOOST_CONCEPT_ASSERT( (bg::concepts::MultiLinestring<multi_t>) );

    multi_t geo;
    bg::read_wkt("MULTILINESTRING((0 0, 5 5, 7 3, 9 5, 10 10),(0 10, 2 8, 4 9))", geo);

    bg::correct(geo);
    BOOST_CHECK_EQUAL(8u, bg::num_points(geo));
    BOOST_CHECK_GT(bg::length(geo), 10.0);

#if defined(BOOST_GEOMETRY_TODO_SUPPORT_CUSTOM_SIMPLIFY)
    // See issue #1133
    multi_t simplified;
    bg::simplify(geo, simplified, 0.1);
    BOOST_CHECK_GT(bg::length(simplified), 10.0);
#endif

    multi_t densified;
    bg::densify(geo, densified, 0.1);
    BOOST_CHECK_GT(bg::num_points(densified), 150u);

    bg::read_wkt("MULTILINESTRING((0 0,  5 5, 5 5, 5 5,   7 3, 9 5, 10 10),(0 10, 2 8, 4 9))", geo);
    BOOST_CHECK_EQUAL(10u, bg::num_points(geo));
    bg::unique(geo);
    BOOST_CHECK_EQUAL(8u, bg::num_points(geo));
}

template <typename P>
void test_two()
{
    using multi_t = cnc_multi_linestring<cnc_linestring<P>>;

    multi_t a;
    a.custom_resize(2);
    fill(a.custom_get(0), {{0, 0}, {5, 5}, {7, 3}, {9, 5}, {10, 10}});
    fill(a.custom_get(1), {{0, 10}, {2, 8}, {4, 9}});

    multi_t b;
    b.custom_resize(2);
    fill(b.custom_get(0), {{0, 0}, {4, 6}, {8, 2}, {9, 8}, {10, 10}});
    fill(b.custom_get(1), {{2, 8}, {4, 9}, {6, 7}});

    BOOST_CHECK_EQUAL(false, bg::crosses(a, b));
    BOOST_CHECK_EQUAL(false, bg::disjoint(a, b));
    BOOST_CHECK_EQUAL(false, bg::equals(a, b));
    BOOST_CHECK_EQUAL(true, bg::intersects(a, b));
    BOOST_CHECK_EQUAL(true, bg::overlaps(a, b));
    BOOST_CHECK_EQUAL(false, bg::touches(a, b));
    BOOST_CHECK_EQUAL(false, bg::relate(a, b, bg::de9im::mask("F0F******")));
    BOOST_CHECK_EQUAL(false, bg::covered_by(a, b));
    BOOST_CHECK_EQUAL(false, bg::within(a, b));

    BOOST_CHECK_LT(bg::comparable_distance(a, b), 1.0);
    BOOST_CHECK_LT(bg::distance(a, b), 1.0);

#if defined(BOOST_GEOMETRY_TODO_SUPPORT_CUSTOM_HAUSDORFF)
    BOOST_CHECK_LT(bg::discrete_hausdorff_distance(a, b), 1.0);
#endif
#if defined(BOOST_GEOMETRY_TODO_SUPPORT_CUSTOM_FRECHET)
    BOOST_CHECK_LT(bg::discrete_frechet_distance(a, b), 1.0);
#endif

#if defined(BOOST_GEOMETRY_TODO_SUPPORT_CUSTOM_INTERSECTION_RESULT)
    using result_t = cnc_multi_linestring<cnc_linestring<P>>;
#else
    using result_t = boost::geometry::model::multi_linestring<boost::geometry::model::linestring<P>>;
#endif

    auto const intersected = linear_intersection<result_t>(a, b);
    BOOST_CHECK_LT(bg::length(intersected), bg::length(a));
    BOOST_CHECK_LT(bg::length(intersected), bg::length(b));

    // There are two linear intersections, and 3 points (degenerated linestrings)
    BOOST_CHECK_EQUAL(5u, bg::num_geometries(intersected));
    BOOST_CHECK_EQUAL(false, bg::is_valid(intersected));

    std::ostringstream svg;
    create_svg(svg, a, b, intersected);

    // "a" as in svg, 10 times and upside down
    std::string svg_a = "0,1000 500,500 700,700 900,500 1000,0";
    BOOST_CHECK_EQUAL(true, svg.str().find(svg_a) != std::string::npos);

    write_svg(svg, "cnc_multi_linestring.svg");
}

template <typename P>
void test_all()
{
    test_const<P>();
    test_mutable<P>();
    test_two<P>();
}

int test_main(int, char* [])
{
    test_all<bg::model::point<double, 2, bg::cs::cartesian> >();

    return 0;
}
