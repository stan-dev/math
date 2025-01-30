# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

foreach(item IN ITEMS
    box_view
    closeable_view
    enumerate_view
    random_access_view
    reversible_closeable
    reversible_view
    segment_view
  )
  boost_geometry_add_unit_test("views" ${item})
endforeach()
