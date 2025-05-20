// Boost.Geometry
// QuickBook Example

// Copyright (c) 2011-2024 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2015 Adam Wulkiewicz, Lodz, Poland.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//[ring
//` Declaration and use of the Boost.Geometry model::ring, modelling the Ring Concept

#include <iostream>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>

namespace bg = boost::geometry;

int main()
{
    using point_t = bg::model::point<double, 2, bg::cs::cartesian>;
    using ring_t = bg::model::ring<point_t>; /*< Default parameters, clockwise, closed ring. >*/

    ring_t ring1; /*< Default-construct a ring. >*/
    ring_t ring2{{0.0, 0.0}, {0.0, 5.0}, {5.0, 5.0}, {5.0, 0.0}, {0.0, 0.0}}; /*< Construct a ring containing four points plus one closing point, using C++11 unified initialization syntax. >*/

    bg::append(ring1, point_t(0.0, 0.0)); /*< Append point. >*/
    bg::append(ring1, point_t(0.0, 5.0));
    bg::append(ring1, point_t(5.0, 5.0));
    bg::append(ring1, point_t(5.0, 0.0));
    bg::append(ring1, point_t(0.0, 0.0));

    double a = bg::area(ring1);

    std::cout << a << std::endl;

    return 0;
}

//]


//[ring_output
/*`
Output:
[pre
25
]
*/
//]
