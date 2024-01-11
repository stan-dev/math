//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/connection.hpp>
#include <boost/mysql/handshake_params.hpp>
#include <boost/mysql/results.hpp>

#include <boost/asio/ssl/verify_mode.hpp>
#include <boost/test/unit_test.hpp>

#include "test_integration/common.hpp"

using namespace boost::mysql::test;
using boost::mysql::common_server_errc;
using boost::mysql::results;

namespace {

auto net_samples_nossl = create_network_samples({
    "tcp_sync_errc",
    "tcp_async_callback",
});

BOOST_AUTO_TEST_SUITE(test_reconnect)

struct reconnect_fixture : network_fixture
{
    void do_query_ok()
    {
        results result;
        conn->execute("SELECT * FROM empty_table", result).get();
        BOOST_TEST(result.rows().empty());
    }
};

BOOST_MYSQL_NETWORK_TEST(reconnect_after_close, reconnect_fixture, net_samples_nossl)
{
    setup(sample.net);

    // Connect and use the connection
    connect();
    do_query_ok();

    // Close
    conn->close().validate_no_error();

    // Reopen and use the connection normally
    connect();
    do_query_ok();
}

BOOST_MYSQL_NETWORK_TEST(reconnect_after_handshake_error, reconnect_fixture, net_samples_nossl)
{
    setup(sample.net);

    // Error during server handshake
    params.set_database("bad_database");
    conn->connect(params).validate_error(
        common_server_errc::er_dbaccess_denied_error,
        {"database", "bad_database"}
    );

    // Reopen with correct parameters and use the connection normally
    params.set_database("boost_mysql_integtests");
    connect();
    do_query_ok();
}

BOOST_AUTO_TEST_SUITE_END()  // test_reconnect

}  // namespace
