// Boost.Geometry (aka GGL, Generic Geometry Library)
// QuickBook Example

// Copyright (c) 2011-2024 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//[assign_points
//` Shows usage of Boost.Geometry's assign, and Boost.Range to assign ranges of a linestring

#include <iostream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>

#include <boost/geometry/geometries/adapted/boost_range/filtered.hpp>

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)

template <typename T>
struct x_between
{
    x_between(T a, T b)
        : fa(a), fb(b)
    {}

    template <typename P>
    bool operator()(P const& p) const
    {
        return boost::geometry::get<0>(p) >= fa
            && boost::geometry::get<0>(p) <= fb;
    }
private :
    T fa, fb;
};


int main()
{
    using ls = boost::geometry::model::linestring<boost::tuple<int, int>>;

    ls line1, line2, line3;

    line1 = {{0, 0}, {2, 3}, {4, 0}, {6, 3}, {8, 0}, {10, 3}, {12, 0}};
    boost::geometry::assign_points(line2, ls({{0, 0}, {2, 2}, {4, 0}, {6, 2}, {8, 0}}));
    boost::geometry::assign_points(line3, line1 | boost::adaptors::filtered(x_between<int>(4, 8))); /*< Boost.Range adaptors can also be used in boost::geometry::assign >*/

    std::cout << "line 1: " << boost::geometry::dsv(line1) << std::endl;
    std::cout << "line 2: " << boost::geometry::dsv(line2) << std::endl;
    std::cout << "line 3: " << boost::geometry::dsv(line3) << std::endl;

    return 0;
}

//]


//[assign_points_output
/*`
Output:
[pre
line 1: ((0, 0), (2, 3), (4, 0), (6, 3), (8, 0), (10, 3), (12, 0))
line 2: ((0, 0), (2, 2), (4, 0), (6, 2), (8, 0))
line 3: ((4, 0), (6, 3), (8, 0))
]
*/
//]
