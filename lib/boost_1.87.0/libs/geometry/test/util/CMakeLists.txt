# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

foreach(item IN ITEMS
    algorithm
    calculation_type
    for_each_coordinate
    math_abs
    math_divide
    math_equals
    math_sqrt
    math_normalize_spheroidal
    promote_integral
    range
    rational
    select_most_precise
    tuples
    write_dsv
  )
  boost_geometry_add_unit_test("util" ${item})
endforeach()
