# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

add_subdirectory(custom_non_copiable)

foreach(item IN ITEMS
    adapted
    boost_array_as_point
    boost_fusion
    # boost_polygon
    boost_range
    boost_tuple
    box
    custom_linestring
    # initialization
    linestring
    multi_linestring
    multi_point
    multi_polygon
    point
    point_xy
    point_xyz
    polygon
    ring
    segment
    infinite_line
  )
  boost_geometry_add_unit_test("geometries" ${item})
endforeach()
