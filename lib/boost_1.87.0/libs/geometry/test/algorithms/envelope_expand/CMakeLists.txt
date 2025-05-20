# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

foreach(item IN ITEMS
    envelope
    envelope_multi
    expand
    expand_on_spheroid
    )
  boost_geometry_add_unit_test("algorithms" ${item})
endforeach()

if (NOT APPLE)
    # The results of these tests vary considerably on Apple/Darwin/arm64 using clang
    foreach(item IN ITEMS
        envelope_on_spheroid
    )
    boost_geometry_add_unit_test("algorithms" ${item})
    endforeach()
endif()