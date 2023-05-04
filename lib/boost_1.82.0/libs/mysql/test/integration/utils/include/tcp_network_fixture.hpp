//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_UTILS_INCLUDE_TCP_NETWORK_FIXTURE_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_UTILS_INCLUDE_TCP_NETWORK_FIXTURE_HPP

#include <boost/mysql/tcp.hpp>

#include <boost/asio/io_context.hpp>

#include "get_endpoint.hpp"
#include "integration_test_common.hpp"
#include "streams.hpp"

namespace boost {
namespace mysql {
namespace test {

struct tcp_network_fixture : network_fixture_base
{
    boost::mysql::tcp_connection conn;

    tcp_network_fixture() : conn(ctx.get_executor()) { conn.set_meta_mode(metadata_mode::full); }

    void connect() { conn.connect(get_endpoint<tcp_socket>(), params); }

    void start_transaction()
    {
        results result;
        conn.query("START TRANSACTION", result);
    }
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif