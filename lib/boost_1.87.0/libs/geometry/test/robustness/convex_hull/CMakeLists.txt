# Boost.Geometry
# Copyright (c) 2024, Oracle and/or its affiliates.
# Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

boost_geometry_add_unit_test("robustness" random_multi_points)

target_include_directories(boost_geometry_robustness_random_multi_points
    PRIVATE
      "${PROJECT_SOURCE_DIR}/test/robustness"
      .)
