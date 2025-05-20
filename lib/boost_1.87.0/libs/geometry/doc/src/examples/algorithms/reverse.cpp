// Boost.Geometry (aka GGL, Generic Geometry Library)
// QuickBook Example

// Copyright (c) 2011-2024 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//[reverse
//` Shows how to reverse a ring or polygon

#include <iostream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)

int main()
{
    using point = boost::tuple<int, int>;
    using polygon = boost::geometry::model::polygon<point>;
    using ring = boost::geometry::model::ring<point>;

    polygon poly;
    boost::geometry::exterior_ring(poly) = {{0, 0}, {0, 9}, {10, 10}, {0, 0}};
    boost::geometry::interior_rings(poly).push_back({{1, 2}, {4, 6}, {2, 8}, {1, 2}});

    double area_before = boost::geometry::area(poly);
    boost::geometry::reverse(poly);
    double area_after = boost::geometry::area(poly);
    std::cout << boost::geometry::dsv(poly) << std::endl;
    std::cout << area_before << " -> " << area_after << std::endl;

    ring r = {{0, 0}, {0, 9}, {8, 8}, {0, 0}};

    area_before = boost::geometry::area(r);
    boost::geometry::reverse(r);
    area_after = boost::geometry::area(r);
    std::cout << boost::geometry::dsv(r) << std::endl;
    std::cout << area_before << " -> " << area_after << std::endl;

    return 0;
}

//]


//[reverse_output
/*`
Output:
[pre
(((0, 0), (10, 10), (0, 9), (0, 0)), ((1, 2), (2, 8), (4, 6), (1, 2)))
38 -> -38
((0, 0), (8, 8), (0, 9), (0, 0))
36 -> -36
]
*/
//]
