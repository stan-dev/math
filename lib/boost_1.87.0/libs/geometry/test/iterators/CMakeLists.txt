# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

foreach(item IN ITEMS
    closing_iterator
    # concatenate_iterator
    ever_circling_iterator
    # flatten_iterator
    point_iterator
    segment_iterator
  )
  boost_geometry_add_unit_test("iterators" ${item})
endforeach()
