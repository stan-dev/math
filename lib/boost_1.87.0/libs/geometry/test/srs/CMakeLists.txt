# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

foreach(item IN ITEMS
    projection
    projection_epsg
    projection_interface_d
    projection_interface_p4
    projection_interface_s
    projection_selftest
    projections
    projections_combined
    projections_static
    srs_transformer
    transformation_epsg
    transformation_interface
  )
  boost_geometry_add_unit_test("srs" ${item})
endforeach()
