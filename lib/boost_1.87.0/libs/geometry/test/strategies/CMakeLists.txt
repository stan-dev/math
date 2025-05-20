# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

foreach(item IN ITEMS
    andoyer
    buffer_join
    buffer_join_geo
    buffer_side_straight
    cross_track
    crossings_multiply
    distance_default_result
    distance_cross_track
    douglas_peucker
    envelope_segment
    franklin
    haversine
    point_in_box
    projected_point
    pythagoras
    pythagoras_point_box
    segment_intersection
    segment_intersection_collinear
    segment_intersection_geo
    segment_intersection_sph
    spherical_side
    side_rounded_input
    thomas
    transform_cs
    transformer
    matrix_transformer
    vincenty
    winding
  )
  boost_geometry_add_unit_test("strategies" ${item})
endforeach()
