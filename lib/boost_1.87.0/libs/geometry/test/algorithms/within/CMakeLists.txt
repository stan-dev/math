# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

# Cartesian
foreach(item IN ITEMS
    within
    within_areal_areal
    within_gc
    within_linear_areal
    within_linear_linear
    within_multi
    within_pointlike_geometry
  )
  boost_geometry_add_unit_test("algorithms" ${item})
endforeach()

# Spherical
foreach(item IN ITEMS
    within_sph
    within_sph_geo
  )
  boost_geometry_add_unit_test("algorithms" ${item})
endforeach()
