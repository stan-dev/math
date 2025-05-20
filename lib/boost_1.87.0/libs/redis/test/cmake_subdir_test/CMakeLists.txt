cmake_minimum_required(VERSION 3.5...3.22)

project(cmake_subdir_test LANGUAGES CXX)

set(BOOST_INCLUDE_LIBRARIES redis)

# Build our dependencies, so the targets Boost::xxx are defined
add_subdirectory(../../../.. boostorg/boost)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE Boost::redis)

include(CTest)
add_test(main main)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C $<CONFIG>)
