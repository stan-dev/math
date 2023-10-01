//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include "test_integration/get_endpoint.hpp"

#include <boost/mysql/connection.hpp>

#include <boost/asio/io_context.hpp>
#include <boost/asio/ip/address_v4.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/test/unit_test.hpp>

#include <cstdlib>

#include "test_integration/safe_getenv.hpp"

namespace {

// Get the endpoint to use for TCP from an environment variable.
// Required as CI MySQL doesn't run on loocalhost
boost::asio::ip::tcp::endpoint get_tcp_valid_endpoint()
{
    std::string hostname = boost::mysql::test::safe_getenv("BOOST_MYSQL_SERVER_HOST", "localhost");
    boost::asio::io_context ctx;
    boost::asio::ip::tcp::resolver resolver(ctx.get_executor());
    auto results = resolver.resolve(hostname, boost::mysql::default_port_string);
    return *results.begin();
}

}  // namespace

boost::asio::ip::tcp::endpoint boost::mysql::test::endpoint_getter<boost::asio::ip::tcp>::operator()()
{
    static auto res = get_tcp_valid_endpoint();
    return res;
}

#ifdef BOOST_ASIO_HAS_LOCAL_SOCKETS
boost::asio::local::stream_protocol::endpoint boost::mysql::test::endpoint_getter<
    boost::asio::local::stream_protocol>::operator()()
{
    return boost::asio::local::stream_protocol::endpoint("/var/run/mysqld/mysqld.sock");
}
#endif
