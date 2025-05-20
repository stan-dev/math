# Copyright 2018, 2019 Peter Dimov
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.5...3.20)

project(cmake_subdir_test LANGUAGES CXX)

add_subdirectory(../.. boostorg/lockfree)

set(deps
    align
    array
    assert
    atomic
    config
    core
    integer
    iterator
    mpl
    parameter
    predef
    static_assert
    tuple
    type_traits
    utility

    # transitive
    bind
    concept_check
    container_hash
    conversion
    describe
    detail
    function
    function_types
    functional
    fusion
    io
    move
    mp11
    optional
    preprocessor
    smart_ptr
    typeof
    throw_exception
    winapi
)

foreach(dep IN LISTS deps)

  add_subdirectory(../../../${dep} boostorg/${dep})

endforeach()

# --target check

add_executable(main ../cmake_install_test/main.cpp)
target_link_libraries(main Boost::lockfree)

enable_testing()
add_test(main main)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C $<CONFIG>)
