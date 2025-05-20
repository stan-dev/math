# Boost.Geometry

# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.

# This file was modified by Oracle on 2024.
# Modifications copyright (c) 2024 Oracle and/or its affiliates.

# Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle

# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

# This takes care of many (but not all) floating point differences on arm64/clang14 on Mac.
# It also lets running the tests much faster.
if (APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffp-contract=fast")
endif()

function(boost_geometry_add_unit_test prefix item)
  set(unit_test_name "boost_geometry_${prefix}_${item}")
  add_executable(${unit_test_name} ${item}.cpp)

  # Add a dependency to Boost.Geometry
  target_link_libraries(${unit_test_name}
    PRIVATE
      Boost::geometry)

  # For unit tests, add a dependency to the unit test framework (in header only mode)
  target_link_libraries(${unit_test_name}
    PRIVATE
      Boost::included_unit_test_framework)

  # Include the main Geometry test folder and the current folder
  target_include_directories(${unit_test_name}
    PRIVATE
      "${PROJECT_SOURCE_DIR}/test"
      "${PROJECT_SOURCE_DIR}/index/test"
      .)

  # To compile with C++14
  target_compile_features(${unit_test_name} PRIVATE cxx_std_14)

  # To be able to run ctest
  set(extra_macro_args ${ARGN})
  if (NOT extra_macro_args STREQUAL "not_run")
    add_test(NAME ${unit_test_name} COMMAND ${unit_test_name})
  endif()

  # Add a dependency to the global tests target
  add_dependencies(tests ${unit_test_name})

  # Inform the caller about the test name. It can then set defines, if necessary.
  set(BOOST_GEOMETRY_UNIT_TEST_NAME ${unit_test_name} PARENT_SCOPE)

endfunction()

if (NOT TARGET tests)
  add_custom_target(tests)
endif()

add_subdirectory(algorithms)
add_subdirectory(arithmetic)
add_subdirectory(concepts)
add_subdirectory(core)
add_subdirectory(cs_undefined)
add_subdirectory(formulas)
add_subdirectory(geometries)
add_subdirectory(io)
add_subdirectory(iterators)
add_subdirectory(policies)
add_subdirectory(robustness)
add_subdirectory(srs)
add_subdirectory(strategies)
add_subdirectory(util)
add_subdirectory(views)
