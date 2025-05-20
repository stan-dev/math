//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_SNIPPETS_GET_ANY_CONNECTION_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_SNIPPETS_GET_ANY_CONNECTION_HPP

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/connect_params.hpp>
#include <boost/mysql/ssl_mode.hpp>

#include <boost/asio/io_context.hpp>

#include "test_common/ci_server.hpp"
#include "test_integration/snippets/credentials.hpp"

namespace boost {
namespace mysql {
namespace test {

struct any_connection_and_context
{
    asio::io_context ctx;
    any_connection conn{ctx.get_executor()};

    any_connection_and_context()
    {
        // Connect
        boost::mysql::connect_params params;
        params.server_address.emplace_host_and_port(get_hostname());
        params.username = mysql_username;
        params.password = mysql_password;
        params.database = "boost_mysql_examples";
        params.ssl = ssl_mode::disable;
        params.multi_queries = true;
        conn.connect(params);
    }
};

inline any_connection& get_any_connection()
{
    static any_connection_and_context obj;
    return obj.conn;
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
