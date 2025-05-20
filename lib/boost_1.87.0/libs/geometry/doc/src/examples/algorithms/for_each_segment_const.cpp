// Boost.Geometry (aka GGL, Generic Geometry Library)
// QuickBook Example

// Copyright (c) 2011-2024 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//[for_each_segment_const
//` Sample using for_each_segment, using a functor to get the minimum and maximum length of a segment in a linestring

#include <iostream>
#include <limits>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

int main()
{
    // Define a type
    using point = boost::geometry::model::d2::point_xy<double>;

    // Create a linestring
    const boost::geometry::model::linestring<point> polyline =
        {{0, 0}, {3, 3}, {5, 1}, {6, 2}, {8, 0}, {4, -4}, {1, -1}, {3, 2}};

    // Iterate over its segments to find the minimum and maximum length
    double min_length = std::numeric_limits<double>::max();
    double max_length = -std::numeric_limits<double>::max();
    
    boost::geometry::for_each_segment(polyline, 
        [&](auto const& s)
        {
            const auto length = static_cast<double>(boost::geometry::length(s));
            min_length = std::min(min_length, length);
            max_length = std::max(max_length, length);
        });

    // Output the results
    std::cout
        << "Min segment length: " << min_length << std::endl
        << "Max segment length: " << max_length << std::endl;

    return 0;
}

//]


//[for_each_segment_const_output
/*`
Output:
[pre
Min segment length: 1.41421
Max segment length: 5.65685
]
*/
//]
