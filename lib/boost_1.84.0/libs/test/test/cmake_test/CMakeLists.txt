# Copyright 2022 Alexander Grund
#
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.5...3.16)

project(cmake_subdir_test LANGUAGES CXX)


# Those 2 should work the same
if(DEFINED BOOST_CI_INSTALL_TEST AND BOOST_CI_INSTALL_TEST)
    message("Using installed Boost")
    find_package(boost_unit_test_framework REQUIRED)
else()
    message("Building Boost")
    add_subdirectory(../.. boostorg/test)
    set(deps
      # Primary dependencies
      algorithm
      assert
      bind
      config
      core
      detail
      exception
      function
      io
      iterator
      mpl
      numeric/conversion
      optional
      preprocessor
      smart_ptr
      static_assert
      type_traits
      utility

      # Secondary dependencies
      array
      concept_check
      range
      regex
      throw_exception
      tuple
      unordered
      typeof
      conversion
      function_types
      fusion
      predef
      move
      container_hash
      functional
      integer
      mp11
      describe
    )

    foreach(dep IN LISTS deps)
      add_subdirectory(../../../${dep} boostorg/${dep})
    endforeach()
endif()

add_executable(main main.cpp)
target_link_libraries(main Boost::unit_test_framework)

enable_testing()
add_test(NAME main COMMAND main)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C $<CONFIG>)
