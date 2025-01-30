# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

foreach(item IN ITEMS
    box
    linestring
    point
    point_xy
    point_xyz
    polygon
    multi_linestring
    multi_point
    multi_polygon
    ring
    segment
  )

  boost_geometry_add_example("geometries" ${item})
  
endforeach()

add_subdirectory(adapted)
add_subdirectory(register)
