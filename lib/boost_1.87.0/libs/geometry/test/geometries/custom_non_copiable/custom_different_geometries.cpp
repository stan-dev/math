// Boost.Geometry
// Unit Test

// Copyright (c) 2023 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// This unit test tests the following
// - if intersection (linear, areal) supports cases where 
//   the input geometry types differ mutually
// - if their point types also differ mutually

#include <boost/geometry.hpp>

#include "helper_functions.hpp"

#include "cnc_container.hpp"
#include "cnc_linestring.hpp"
#include "cnc_ring.hpp"
#include "cnc_polygon.hpp"
#include "cnc_multi_linestring.hpp"
#include "cnc_multi_polygon.hpp"

#include "adapt_cnc_container.hpp"
#include "adapt_cnc_linestring.hpp"
#include "adapt_cnc_ring.hpp"
#include "adapt_cnc_polygon.hpp"
#include "adapt_cnc_multi_linestring.hpp"
#include "adapt_cnc_multi_polygon.hpp"

#include <geometry_test_common.hpp>

#include <sstream>

namespace bg = boost::geometry;

template <typename P1, typename P2>
void test_multi_poly_different_types()
{
    using multi_1_t = cnc_multi_polygon<cnc_polygon<P1>>;
    using multi_2_t = boost::geometry::model::multi_polygon<boost::geometry::model::polygon<P2>>;

    multi_1_t a;
    a.custom_resize(1);
    a.custom_get(0).custom_int().custom_resize(2);
    fill(a.custom_get(0).custom_ext(), {{0, 0}, {0, 5}, {5, 5}, {5, 0}, {0, 0}});
    fill(a.custom_get(0).custom_int().custom_get(0), {{1, 1}, {2, 1}, {2, 2}, {1, 2}, {1, 1}});
    fill(a.custom_get(0).custom_int().custom_get(1), {{3, 3}, {4, 3}, {4, 4}, {3, 4}, {3, 3}});

    multi_2_t b;
    bg::read_wkt("MULTIPOLYGON(((1 1,1 2,3 2,3 1,1 1)),((3 3,3 4,5 4,5 3,3 3)))", b);

    auto const intersected = areal_intersection<multi_2_t>(a, b);
    BOOST_CHECK_LT(bg::area(intersected), bg::area(a));
    BOOST_CHECK_LT(bg::area(intersected), bg::area(b));

    std::ostringstream svg;
    create_svg(svg, a, b, intersected);
    write_svg(svg, "cnc_multi_polygon_different.svg");
}

template <typename P1, typename P2>
void test_multi_linestring_different_types()
{
    using multi_1_t = cnc_multi_linestring<cnc_linestring<P1>>;
    using multi_2_t = boost::geometry::model::multi_linestring<boost::geometry::model::linestring<P2>>;

    multi_1_t a;
    multi_2_t b;
    bg::read_wkt("MULTILINESTRING((0 0, 5 5, 7 3, 9 5, 10 10),(0 10, 2 8, 4 9))", a);
    bg::read_wkt("MULTILINESTRING((0 0, 4 6, 8 2, 9 8, 10 10), (2 8, 4 9, 6 7))", b);

    auto const intersected = linear_intersection<multi_2_t>(a, b);
    BOOST_CHECK_LT(bg::length(intersected), bg::length(a));
    BOOST_CHECK_LT(bg::length(intersected), bg::length(b));

    std::ostringstream svg;
    create_svg(svg, a, b, intersected);
    write_svg(svg, "cnc_multi_linestring_different.svg");
}

template <typename P1, typename P2>
void test_all()
{
    test_multi_poly_different_types<P1, P2>();
    test_multi_linestring_different_types<P1, P2>();
}

int test_main(int, char* [])
{
    test_all
        <
            bg::model::point<double, 2, bg::cs::cartesian>,
            bg::model::point<float, 2, bg::cs::cartesian>
        >();

    return 0;
}
