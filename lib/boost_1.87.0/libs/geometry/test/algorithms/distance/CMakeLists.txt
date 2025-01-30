# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

foreach(item IN ITEMS
    distance
    distance_ca_ar_ar
    distance_ca_l_ar
    distance_ca_l_l
    distance_ca_pl_ar
    distance_ca_pl_l
    distance_ca_pl_pl
    distance_se_geo_ar_ar
    distance_se_geo_l_ar
    distance_se_geo_l_l
    distance_se_geo_pl_ar
    distance_geo_pl_l
    distance_se_geo_pl_pl
    distance_se_pl_l
    distance_se_pl_pl
  )
  boost_geometry_add_unit_test("algorithms" ${item})
endforeach()
