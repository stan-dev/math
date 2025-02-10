#
# Copyright (c) 2023 Dmitry Arkhipov (grisumbras@yandex.ru)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/json
#

cmake_minimum_required(VERSION 3.8...3.16)

project(myproject LANGUAGES CXX)

add_subdirectory(../../../align           boostorg/align)
add_subdirectory(../../../assert          boostorg/assert)
add_subdirectory(../../../config          boostorg/config)
add_subdirectory(../../../container       boostorg/container)
add_subdirectory(../../../container_hash  boostorg/container_hash)
add_subdirectory(../../../core            boostorg/core)
add_subdirectory(../../../describe        boostorg/describe)
add_subdirectory(../../../endian          boostorg/endian)
add_subdirectory(../../../intrusive       boostorg/intrusive)
add_subdirectory(../../../json            boostorg/json)
add_subdirectory(../../../move            boostorg/move)
add_subdirectory(../../../mp11            boostorg/mp11)
add_subdirectory(../../../predef          boostorg/predef)
add_subdirectory(../../../system          boostorg/system)
add_subdirectory(../../../throw_exception boostorg/throw_exception)
add_subdirectory(../../../static_assert   boostorg/static_assert)
add_subdirectory(../../../variant2        boostorg/variant2)
add_subdirectory(../../../winapi          boostorg/winapi)

add_executable(mytest main.cpp)
target_link_libraries(mytest PRIVATE Boost::json)

enable_testing()
add_test(mytest mytest)
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure --no-tests=error -C $<CONFIG>)
add_dependencies(check mytest)
