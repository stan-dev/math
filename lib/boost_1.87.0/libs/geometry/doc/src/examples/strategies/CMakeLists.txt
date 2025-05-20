# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

foreach(item IN ITEMS
    buffer_point_circle
    buffer_geographic_point_circle
    buffer_point_square
    buffer_join_round
    buffer_join_miter
    buffer_end_round
    buffer_end_flat
    buffer_distance_symmetric
    buffer_distance_asymmetric
  )

  boost_geometry_add_example("strategies" ${item})
  
endforeach()

