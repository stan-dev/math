cmake_minimum_required(VERSION 3.5...3.22)

project(cmake_install_test LANGUAGES CXX)

find_package(boost_redis REQUIRED)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE Boost::redis)

include(CTest)
add_test(main main)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C $<CONFIG>)
