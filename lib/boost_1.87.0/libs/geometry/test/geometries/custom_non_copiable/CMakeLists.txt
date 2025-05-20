# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

foreach(item IN ITEMS
    custom_ring
    custom_polygon
    custom_multi_polygon
    custom_linestring
    custom_multi_linestring
    custom_different_geometries
  )
  boost_geometry_add_unit_test("geometries_cnc" ${item})
endforeach()
