// Boost.Geometry (aka GGL, Generic Geometry Library)
// QuickBook Example

// Copyright (c) 2011-2024 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//[assign_point_to_index
//` Shows how to assign the lower-left or upper-right point from a box

#include <iostream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

using namespace boost::geometry;

int main()
{
    using point = model::d2::point_xy<int>;
    using box = model::box<point>;

    point lower_left(0, 0), upper_right(2, 2);

    box b;
    detail::assign_point_to_index<0>(lower_left, b);
    detail::assign_point_to_index<1>(upper_right, b);
    std::cout << "box: " << dsv(b) << std::endl;

    return 0;
}

//]


//[assign_point_to_index_output
/*`
Output:
[pre
box: ((0, 0), (2, 2))
]
*/
//]
