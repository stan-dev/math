# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

foreach(item IN ITEMS
          union
          union_multi
          union_other_types
        )
  boost_geometry_add_unit_test("algorithms" ${item})
  target_compile_definitions(${BOOST_GEOMETRY_UNIT_TEST_NAME} PRIVATE BOOST_GEOMETRY_TEST_ONLY_ONE_TYPE)
endforeach()

foreach(item IN ITEMS
          union_aa_geo
          union_aa_sph
          union_gc
          union_linear_linear
          union_pl_pl
          union_issues
          union_tupled
        )
  boost_geometry_add_unit_test("algorithms" ${item})
endforeach()
