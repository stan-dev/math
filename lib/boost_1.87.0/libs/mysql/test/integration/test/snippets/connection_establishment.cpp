//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/connect_params.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/io_context_fixture.hpp"

namespace mysql = boost::mysql;
namespace asio = boost::asio;
using namespace boost::mysql::test;

namespace {

BOOST_FIXTURE_TEST_CASE(section_connection_establishment, io_context_fixture)
{
    {
        //[section_connection_establishment_max_buffer_size
        // Increase the max buffer size to 512MB.
        // This allows reading individual rows as big as 512MB.
        // This is only required if each individual row is extremely big,
        // and is not required for many smaller rows.
        mysql::any_connection_params conn_params;
        conn_params.max_buffer_size = 0x20000000;

        // Create the connection
        mysql::any_connection conn(ctx, conn_params);

        // Connect and use the connection normally
        //]
    }
}

}  // namespace
