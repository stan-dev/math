//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_SNIPPETS_GET_CONNECTION_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_SNIPPETS_GET_CONNECTION_HPP

#include <boost/mysql/tcp.hpp>

#include <boost/asio/io_context.hpp>

#include "test_common/ci_server.hpp"
#include "test_integration/snippets/credentials.hpp"

namespace boost {
namespace mysql {
namespace test {

struct connection_and_context
{
    asio::io_context ctx;
    tcp_connection conn{ctx.get_executor()};

    connection_and_context()
    {
        // Resolve the hostname to get a collection of endpoints
        boost::asio::ip::tcp::resolver resolver(ctx.get_executor());
        auto endpoints = resolver.resolve(get_hostname(), boost::mysql::default_port_string);

        // Connect
        boost::mysql::handshake_params params(mysql_username, mysql_password, "boost_mysql_examples");
        conn.connect(*endpoints.begin(), params);
    }
};

inline tcp_connection& get_connection()
{
    static connection_and_context obj;
    return obj.conn;
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
