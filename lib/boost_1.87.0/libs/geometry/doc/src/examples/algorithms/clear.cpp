// Boost.Geometry (aka GGL, Generic Geometry Library)
// QuickBook Example

// Copyright (c) 2011-2024 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//[clear
//` Shows how to clear a ring or polygon

#include <iostream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)

int main()
{
    using point = boost::tuple<float, float>;
    using polygon = boost::geometry::model::polygon<point>;
    using ring = boost::geometry::model::ring<point>;

    // Create a square as exterior ring, with a triangle as interior ring
    polygon poly{{{0, 0}, {0, 10}, {10, 10}, {10, 0}, {0, 0}}, {{{1, 2}, {8, 2}, {4, 6}, {1, 2}}}};

    std::cout << boost::geometry::dsv(poly) << std::endl;
    boost::geometry::clear(poly);
    std::cout << boost::geometry::dsv(poly) << std::endl;

    // Create a triangle
    ring r{{0, 0}, {0, 9}, {8, 8}, {0, 0}};

    std::cout << boost::geometry::dsv(r) << std::endl;
    boost::geometry::clear(r);
    std::cout << boost::geometry::dsv(r) << std::endl;

    return 0;
}

//]


//[clear_output
/*`
Output:
[pre
(((0, 0), (0, 10), (10, 10), (10, 0), (0, 0)), ((1, 2), (8, 2), (4, 6), (1, 2)))
(())
((0, 0), (0, 9), (8, 8), (0, 0))
()
]
*/
//]
