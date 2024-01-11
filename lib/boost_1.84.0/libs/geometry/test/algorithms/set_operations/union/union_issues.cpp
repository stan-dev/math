// Boost.Geometry

// Unit Test

// Copyright (c) 2023 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// This unit test is verifying some issues which do not use WKT,
// or are not expressable as WKT (for example by the use of an epsilon)

#include <geometry_test_common.hpp>

#include <boost/geometry/algorithms/area.hpp>
#include <boost/geometry/algorithms/union.hpp>
#include <boost/geometry/geometries/geometries.hpp>

#include <limits>

template <typename T>
void issue_1103()
{
    using point_t = bg::model::d2::point_xy<T>;
    using polygon_t = bg::model::polygon<point_t, true, true>;

    polygon_t poly1;
    bg::append(poly1, point_t(1, 1));
    bg::append(poly1, point_t(1, std::numeric_limits<T>::epsilon()));
    bg::append(poly1, point_t(0, std::numeric_limits<T>::epsilon()));
    bg::append(poly1, point_t(0, 1));
    bg::append(poly1, point_t(1, 1));

    polygon_t poly2;
    bg::append(poly2, point_t(1, 0));
    bg::append(poly2, point_t(1, -1));
    bg::append(poly2, point_t(0, -1));
    bg::append(poly2, point_t(0, 0));
    bg::append(poly2, point_t(1, 0));

    bg::model::multi_polygon<polygon_t> result;
    bg::union_(poly1, poly2, result);

    // Verify result. Before commit b1bebca the result was empty.
    BOOST_CHECK_EQUAL(1, static_cast<int>(boost::size(result)));
    BOOST_CHECK_CLOSE(2.0, bg::area(result), 0.0001);
}

int test_main(int, char* [])
{
    issue_1103<double>();
    issue_1103<long double>();

    return 0;
}
