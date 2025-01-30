# Boost.Geometry
# Copyright (c) 2024, Oracle and/or its affiliates.
# Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

foreach(item IN ITEMS
    rtree_contains_point
    rtree_epsilon
    rtree_insert_remove
    rtree_intersects_geom
    rtree_move_pack
    rtree_non_cartesian
    rtree_values
    #compile-fail rtree_values_invalid
  )
  boost_geometry_add_unit_test("index" ${item})
endforeach()
