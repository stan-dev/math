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
#include "cnc_multi_polygon.hpp"

#include "adapt_cnc_container.hpp"
#include "adapt_cnc_ring.hpp"
#include "adapt_cnc_polygon.hpp"

#include <geometry_test_common.hpp>

#include <sstream>

namespace bg = boost::geometry;

template <typename P>
void test_const()
{
    using polygon_t = cnc_polygon<P>;
    using point_t = typename bg::point_type<polygon_t>::type;

    boost::ignore_unused<point_t>();

    BOOST_CONCEPT_ASSERT( (bg::concepts::ConstPolygon<polygon_t>) );

    polygon_t geo;
    geo.custom_int().custom_resize(2);
    fill(geo.custom_ext(), {{0, 0}, {0, 5}, {5, 5}, {5, 0}, {0, 0}});
    fill(geo.custom_int().custom_get(0), {{1, 1}, {2, 1}, {2, 2}, {1, 2}, {1, 1}});
    fill(geo.custom_int().custom_get(1), {{3, 3}, {4, 3}, {4, 4}, {3, 4}, {3, 3}});

#if defined(TEST_FAIL_CNC)
    // This should NOT work.
    auto copy = geo;
#endif

    BOOST_CHECK_EQUAL(15u, bg::num_points(geo));
    BOOST_CHECK_EQUAL(false, bg::is_empty(geo));
    BOOST_CHECK_EQUAL(false, bg::is_convex(geo));
    BOOST_CHECK_EQUAL(true, bg::is_simple(geo));
    BOOST_CHECK_EQUAL(true, bg::is_valid(geo));
    BOOST_CHECK_EQUAL(12u, bg::num_segments(geo));
    BOOST_CHECK_EQUAL(2u, bg::num_interior_rings(geo));

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
    out << "POLYGON: " << bg::wkt(geo)
        << " centroid: " << bg::wkt(centroid(geo))
        << " envelope: " << bg::wkt(envelope(geo))
        << " hull: " << bg::wkt(hull(geo))
        << " surface: " << bg::wkt(point_on_surface(geo));

    // Test Boost.Range behavior over interior rings
    auto const& ring = bg::exterior_ring(geo);
    out << "RINGS: " << bg::wkt(ring) << bg::dsv(ring);

    // This should NOT compile
#if defined(TEST_FAIL_CNC)
    auto rings = bg::interior_rings(geo);
#endif

    auto const& rings = bg::interior_rings(geo);
    for (auto it = boost::begin(rings); it != boost::end(rings); ++it)
    {
        out << " " << bg::wkt(*it) << " " << bg::dsv(*it);
    }
    out << std::endl;

#if defined(TEST_FAIL_CNC)
    // Test std::begin/end (does not work)
    for (auto it = std::begin(rings); it != std::end(rings); ++it)
    {
        out << " " << bg::wkt(*it) << bg::dsv(*it);
    }
    out << std::endl;

    // Test range based for loops (does not work)
    for (auto const& ring : rings)
    {
        out << " " << bg::wkt(ring) << " " << bg::dsv(ring);
    }
    out << std::endl;
#endif
    auto const test_string = out.str();
    auto const test_length = test_string.length();
    BOOST_CHECK_MESSAGE(test_length > 300,
        "detected: '" << test_string << "' (" << test_length << ")");
}

template <typename P>
void test_mutable()
{
    using polygon_t = cnc_polygon<P>;
    using point_t = typename bg::point_type<polygon_t>::type;

    boost::ignore_unused<point_t>();

    BOOST_CONCEPT_ASSERT( (bg::concepts::Polygon<polygon_t>) );

    polygon_t geo;
    // TODO: support WKT reading for polygons with non-copyable rings
    // For now read it into the exterior ring
    bg::read_wkt("POLYGON((0 0,0 4,4 4,4 0))", geo.custom_ext());

    bg::correct(geo);
    BOOST_CHECK_EQUAL(5u, bg::num_points(geo));
    BOOST_CHECK_GT(bg::area(geo), 10.0);

#if defined(BOOST_GEOMETRY_TODO_SUPPORT_CUSTOM_SIMPLIFY)
    // See issue #1133
    polygon_t simplified;
    bg::simplify(geo, simplified, 0.1);
    BOOST_CHECK_GT(bg::area(simplified), 10.0);
#endif

    polygon_t densified;
    bg::densify(geo, densified, 0.1);
    BOOST_CHECK_GT(bg::num_points(densified), 150u);

    bg::read_wkt("POLYGON((0 0,0 4,4 4,4 4,4 4,4 0,0 0))", geo.custom_ext());
    bg::unique(geo);
    BOOST_CHECK_EQUAL(5u, bg::num_points(geo));

    auto const converted = convert_to<cnc_ring<P>>(geo);
    BOOST_CHECK_EQUAL(5u, bg::num_points(converted));
}

template <typename P>
void test_two()
{
    using polygon_t = cnc_polygon<P>;

    polygon_t a;
    a.custom_int().custom_resize(2);
    fill(a.custom_ext(), {{0, 0}, {0, 5}, {5, 5}, {5, 0}, {0, 0}});
    fill(a.custom_int().custom_get(0), {{1, 1}, {2, 1}, {2, 2}, {1, 2}, {1, 1}});
    fill(a.custom_int().custom_get(1), {{3, 3}, {4, 3}, {4, 4}, {3, 4}, {3, 3}});
    polygon_t b;
    fill(b.custom_ext(), {{1, 1}, {1, 5}, {5, 5}, {5, 1}, {1, 1}});

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
    using result_t = cnc_multi_polygon<cnc_polygon<P>>;
#else
    using result_t = boost::geometry::model::multi_polygon<boost::geometry::model::polygon<P>>;
#endif

    auto const intersected = areal_intersection<result_t>(a, b);
    BOOST_CHECK_LT(bg::area(intersected), bg::area(a));
    BOOST_CHECK_LT(bg::area(intersected), bg::area(b));

    std::ostringstream svg;
    create_svg(svg, a, b, intersected);

    // "b" as in svg, 10 times and upside down
    std::string svg_b = "M 200,800 L 200,0 L 1000,0 L 1000,800 L 200,800 z";
    BOOST_CHECK_EQUAL(true, svg.str().find(svg_b) != std::string::npos);

    write_svg(svg, "cnc_polygon.svg");
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
