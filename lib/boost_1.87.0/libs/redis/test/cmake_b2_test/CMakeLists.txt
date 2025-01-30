cmake_minimum_required(VERSION 3.5...3.22)

project(boost_redis_b2_test LANGUAGES CXX)

find_package(Boost REQUIRED COMPONENTS headers)
find_package(Threads REQUIRED)
find_package(OpenSSL REQUIRED)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE Boost::headers Threads::Threads OpenSSL::Crypto OpenSSL::SSL)

include(CTest)
add_test(NAME main COMMAND main)
