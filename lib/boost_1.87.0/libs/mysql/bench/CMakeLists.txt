#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#


add_executable(
    boost_mysql_bench_connection_pool
    connection_pool.cpp
)

target_link_libraries(
    boost_mysql_bench_connection_pool
    PUBLIC
    boost_mysql_compiled
)

boost_mysql_common_target_settings(boost_mysql_bench_connection_pool)