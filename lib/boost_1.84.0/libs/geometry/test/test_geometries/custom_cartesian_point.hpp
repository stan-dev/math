// Boost.Geometry

// Copyright (c) 2022 Oracle and/or its affiliates.

// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/geometry/geometries/geometries.hpp>

struct custom_cartesian_point {
    double x, y;
    custom_cartesian_point(double x = 0, double y = 0) : x(x), y(y) {}
};

BOOST_GEOMETRY_REGISTER_POINT_2D(custom_cartesian_point, double,
                                 boost::geometry::cs::cartesian, x, y);

using custom_cartesian_line = boost::geometry::model::linestring<custom_cartesian_point>;
using custom_cartesian_polygon = boost::geometry::model::polygon<custom_cartesian_point>;