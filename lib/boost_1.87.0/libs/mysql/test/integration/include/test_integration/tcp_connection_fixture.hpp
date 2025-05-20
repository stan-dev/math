//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_TCP_CONNECTION_FIXTURE_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_TCP_CONNECTION_FIXTURE_HPP

#include <boost/mysql/handshake_params.hpp>
#include <boost/mysql/tcp.hpp>

#include <boost/assert/source_location.hpp>

#include "test_common/io_context_fixture.hpp"
#include "test_common/source_location.hpp"

namespace boost {
namespace mysql {
namespace test {

struct tcp_connection_fixture : io_context_fixture
{
    tcp_connection conn;

    tcp_connection_fixture();
    ~tcp_connection_fixture();

    void connect(source_location loc = BOOST_MYSQL_CURRENT_LOCATION);
    void connect(const handshake_params&, source_location loc = BOOST_MYSQL_CURRENT_LOCATION);
};

asio::ip::tcp::endpoint get_tcp_endpoint();

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif