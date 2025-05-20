// Boost.Geometry
// Unit Test

// Copyright (c) 2023 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/geometry.hpp>

#include "helper_functions.hpp"

#include "cnc_container.hpp"
#include "cnc_ring.hpp"
#include "cnc_polygon.hpp"

#include "adapt_cnc_container.hpp"
#include "adapt_cnc_ring.hpp"
#include "adapt_cnc_polygon.hpp"

#include <geometry_test_common.hpp>

#include <sstream>

namespace bg = boost::geometry;

template <typename P>
void test_const()
{
    using ring_t = cnc_ring<P>;
    using point_t = typename bg::point_type<ring_t>::type;

    boost::ignore_unused<point_t>();

    BOOST_CONCEPT_ASSERT( (bg::concepts::ConstRing<ring_t>) );

    ring_t geo;
    fill(geo, {{0, 0}, {0, 4}, {4, 4}, {4, 0}, {0, 0}});

#if defined(TEST_FAIL_CNC)
    // This should NOT work.
    auto copy = geo;
#endif

    BOOST_CHECK_EQUAL(5u, boost::size(geo));
    BOOST_CHECK_EQUAL(5u, bg::num_points(geo));
    BOOST_CHECK_EQUAL(false, bg::is_empty(geo));
    BOOST_CHECK_EQUAL(true, bg::is_convex(geo));
    BOOST_CHECK_EQUAL(true, bg::is_simple(geo));
    BOOST_CHECK_EQUAL(true, bg::is_valid(geo));
    BOOST_CHECK_EQUAL(4u, bg::num_segments(geo));
    BOOST_CHECK_EQUAL(0u, bg::num_interior_rings(geo));

    // Check floating point properties, they have to compile, matching exactly is done elsewhere
    BOOST_CHECK_GT(bg::area(geo), 10.0);
    BOOST_CHECK_GT(bg::perimeter(geo), 10.0);

#if defined(BOOST_GEOMETRY_TODO_SUPPORT_CUSTOM_BUFFER_RESULT)
    using result_t = cnc_multi_polygon<cnc_polygon<P>>;
#else
    using result_t = boost::geometry::model::multi_polygon<boost::geometry::model::polygon<P>>;
#endif
    auto const buffered = buffer<result_t>(geo, 0.1);
    BOOST_CHECK_GT(bg::area(buffered), bg::area(geo));
    BOOST_CHECK_EQUAL(true, bg::within(geo, buffered));

    std::ostringstream out;
    out << "RING: " << bg::wkt(geo)
        << " centroid: " << bg::wkt(centroid(geo))
        << " envelope: " << bg::wkt(envelope(geo))
        << " hull: " << bg::wkt(hull(geo))
        << " surface: " << bg::wkt(point_on_surface(geo));

    // Test Boost.Range behavior
    for (auto it = boost::begin(geo); it != boost::end(geo); ++it)
    {
        out << " " << bg::wkt(*it);
    }
    out << std::endl;

#if 0
    // Test std::begin/end (does not work)
    for (auto it = std::begin(geo); it != std::end(geo); ++it)
    {
        out << " " << bg::wkt(*it);
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
    using ring_t = cnc_ring<P>;
    using point_t = typename bg::point_type<ring_t>::type;

    boost::ignore_unused<point_t>();

    BOOST_CONCEPT_ASSERT( (bg::concepts::Ring<ring_t>) );

    ring_t geo;
    bg::read_wkt("POLYGON((0 0,0 4,4 4,4 0))", geo);

    bg::correct(geo);
    BOOST_CHECK_EQUAL(5u, bg::num_points(geo));
    BOOST_CHECK_GT(bg::area(geo), 10.0);

#if defined(BOOST_GEOMETRY_TODO_SUPPORT_CUSTOM_SIMPLIFY)
    // See issue #1133
    ring_t simplified;
    bg::simplify(geo, simplified, 0.1);
    BOOST_CHECK_GT(bg::area(simplified), 10.0);
#endif

    ring_t densified;
    bg::densify(geo, densified, 0.1);
    BOOST_CHECK_GT(bg::num_points(densified), 150u);

    bg::read_wkt("POLYGON((0 0,0 4,4 4,4 4,4 4,4 0,0 0))", geo);
    bg::unique(geo);
    BOOST_CHECK_EQUAL(5u, bg::num_points(geo));

    auto const converted = convert_to<cnc_polygon<P>>(geo);
    BOOST_CHECK_EQUAL(5u, bg::num_points(converted));
}

template <typename P>
void test_two()
{
    using ring_t = cnc_ring<P>;

    ring_t a;
    fill(a, {{0, 0}, {0, 4}, {4, 4}, {4, 0}, {0, 0}});

    ring_t b;
    fill(b, {{1, 1}, {1, 5}, {5, 5}, {5, 1}, {1, 1}});

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
    // See issue #1134
    BOOST_CHECK_LT(bg::discrete_hausdorff_distance(a, b), 1.0);
#endif
#if defined(BOOST_GEOMETRY_TODO_SUPPORT_CUSTOM_FRECHET)
    // See issue #1134
    BOOST_CHECK_LT(bg::discrete_frechet_distance(a, b), 1.0);
#endif

#if defined(BOOST_GEOMETRY_TODO_SUPPORT_CUSTOM_INTERSECTION_RESULT)
    # See issue #1135
    using result_t = cnc_multi_polygon<cnc_polygon<P>>;
#else
    using result_t = boost::geometry::model::multi_polygon<boost::geometry::model::polygon<P>>;
#endif

    auto const intersected = areal_intersection<result_t>(a, b);
    BOOST_CHECK_LT(bg::area(intersected), bg::area(a));
    BOOST_CHECK_LT(bg::area(intersected), bg::area(b));
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
