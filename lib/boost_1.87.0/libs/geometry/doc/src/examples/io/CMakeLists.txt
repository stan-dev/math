# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

foreach(item IN ITEMS
    svg
    svg_mapper_custom
    wkt
    read_wkt
    to_wkt
    from_wkt
  )

  boost_geometry_add_example("io" ${item})
  
endforeach()

