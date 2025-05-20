# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

foreach(item IN ITEMS
    buffer
    buffer_gc
    buffer_variable_width
  )
  boost_geometry_add_unit_test("algorithms" ${item})
endforeach()

# Tests for all geometry types, cartesian, for one floating point type
foreach(item IN ITEMS
    buffer_point
    buffer_linestring
    buffer_ring
    buffer_polygon
    buffer_multi_point
    buffer_multi_polygon
    buffer_multi_linestring
  )
  boost_geometry_add_unit_test("algorithms" ${item})
  target_compile_definitions(${BOOST_GEOMETRY_UNIT_TEST_NAME} PRIVATE BOOST_GEOMETRY_TEST_ONLY_ONE_TYPE)
endforeach()

# Tests for geographic and spherical
foreach(item IN ITEMS
    buffer_point_geo
    buffer_polygon_geo
    buffer_multi_linestring_geo
    buffer_multi_polygon_geo
    buffer_geo_spheroid
  )
  boost_geometry_add_unit_test("algorithms" ${item})
  target_compile_definitions(${BOOST_GEOMETRY_UNIT_TEST_NAME} PRIVATE BOOST_GEOMETRY_TEST_ONLY_ONE_TYPE)
endforeach()

# Other tests
foreach(item IN ITEMS
    buffer_with_strategies
    buffer_piece_border
    buffer_linestring_aimes
  )
  boost_geometry_add_unit_test("algorithms" ${item})
  target_compile_definitions(${BOOST_GEOMETRY_UNIT_TEST_NAME} PRIVATE BOOST_GEOMETRY_TEST_ONLY_ONE_TYPE)
endforeach()
