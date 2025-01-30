# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

# Cartesian
foreach(item IN ITEMS
    assemble
    copy_segment_point
    get_clusters
    get_distance_measure
    get_ring
    get_turn_info
    get_turns
    get_turns_const
    get_turns_areal_areal
    get_turns_linear_areal
    get_turns_linear_linear
    overlay
    sort_by_side_basic
    sort_by_side
    relative_order
    select_rings
    self_intersection_points
  )
  boost_geometry_add_unit_test("algorithms" ${item})
endforeach()

# Spherical
foreach(item IN ITEMS
    get_turns_areal_areal_sph
    get_turns_linear_areal_sph
    get_turns_linear_linear_sph
    get_turns_linear_linear_geo
  )
  boost_geometry_add_unit_test("algorithms" ${item})
endforeach()
