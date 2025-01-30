# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

# Cartesian
foreach(item IN ITEMS
    relate_const_custom
    relate_areal_areal
    relate_linear_areal
    relate_linear_linear
    relate_pointlike_geometry
  )
  boost_geometry_add_unit_test("algorithms" ${item})
endforeach()

# Spherical
foreach(item IN ITEMS
    relate_areal_areal_sph
    relate_linear_areal_sph
    relate_linear_linear_sph
  )
  boost_geometry_add_unit_test("algorithms" ${item})
endforeach()
