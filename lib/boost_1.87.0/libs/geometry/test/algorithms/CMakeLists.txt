# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

add_subdirectory(area)
add_subdirectory(buffer)
add_subdirectory(covered_by)
add_subdirectory(convex_hull)
add_subdirectory(closest_points)
add_subdirectory(crosses)
add_subdirectory(detail)
add_subdirectory(disjoint)
add_subdirectory(distance)
add_subdirectory(envelope_expand)
add_subdirectory(equals)
add_subdirectory(intersects)
add_subdirectory(length)
add_subdirectory(overlaps)
add_subdirectory(overlay)
add_subdirectory(perimeter)
add_subdirectory(relate)
add_subdirectory(set_operations)
add_subdirectory(similarity)
add_subdirectory(touches)
add_subdirectory(within)

foreach(item IN ITEMS
    append
    assign
    azimuth
    centroid
    centroid_multi
    comparable_distance
    convert
    convert_multi
    correct
    correct_multi
    correct_closure
    densify
    for_each
    for_each_multi
    is_convex
    is_empty
    is_simple
    is_simple_geo
    is_valid
    is_valid_failure
    is_valid_geo
    line_interpolate
    make
    maximum_gap
    merge_elements
    num_geometries
    num_geometries_multi
    num_interior_rings
    num_interior_rings_multi
    num_points
    num_points_multi
    num_segments
    point_on_surface
    remove_spikes
    reverse
    reverse_multi
    simplify
    simplify_multi
    transform
    transform_multi
    unique
    unique_multi
  )
  boost_geometry_add_unit_test("algorithms" ${item})
endforeach()
