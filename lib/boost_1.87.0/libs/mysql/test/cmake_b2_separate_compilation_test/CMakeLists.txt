#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# This is included in 02_integrating.qbk. If you change it, update these files, too.

cmake_minimum_required(VERSION 3.5...3.22)

project(boost_mysql_example LANGUAGES CXX)

find_package(Boost REQUIRED COMPONENTS headers charconv)
find_package(Threads REQUIRED)
find_package(OpenSSL REQUIRED)

add_executable(
    main
    # Contains Boost.MySQL sources via #include <boost/mysql/src.hpp>
    boost_mysql.cpp
    # List any other .cpp your exe has here
    main.cpp
)
target_link_libraries(main PRIVATE Boost::charconv Threads::Threads OpenSSL::Crypto OpenSSL::SSL)

# We need to define BOOST_MYSQL_SEPARATE_COMPILATION in any code using Boost.MySQL in separate-build mode
target_compile_definitions(main PRIVATE BOOST_MYSQL_SEPARATE_COMPILATION)

include(CTest)
add_test(NAME main COMMAND main)
