//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/connect_params.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/ssl_mode.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/asio/bind_cancellation_slot.hpp>
#include <boost/asio/bind_executor.hpp>
#include <boost/asio/cancellation_signal.hpp>
#include <boost/asio/cancellation_type.hpp>
#include <boost/asio/error.hpp>
#include <boost/asio/post.hpp>
#include <boost/core/span.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "test_common/create_basic.hpp"
#include "test_common/network_result.hpp"
#include "test_integration/any_connection_fixture.hpp"
#include "test_integration/connect_params_builder.hpp"
#include "test_integration/server_features.hpp"
#include "test_integration/spotchecks_helpers.hpp"
#include "test_integration/tcp_connection_fixture.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;
namespace asio = boost::asio;
using boost::test_tools::per_element;

namespace {

BOOST_AUTO_TEST_SUITE(test_reconnect)

auto connection_samples = network_functions_connection::sync_and_async();
auto any_samples = network_functions_any::sync_and_async();
auto any_samples_grid = boost::unit_test::data::make(any_samples) *
                        boost::unit_test::data::make({ssl_mode::disable, ssl_mode::require});

// Old connection can reconnect after close if the stream is not SSL
BOOST_DATA_TEST_CASE(reconnect_after_close_connection, connection_samples)
{
    netfn_fixture_connection fix(sample);

    // Connect and use the connection
    results r;
    fix.net.connect(fix.conn, get_tcp_endpoint(), connect_params_builder().build_hparams())
        .validate_no_error();
    fix.net.execute_query(fix.conn, "SELECT * FROM empty_table", r).validate_no_error();

    // Close
    fix.net.close(fix.conn).validate_no_error();

    // Reopen and use the connection normally
    fix.net.connect(fix.conn, get_tcp_endpoint(), connect_params_builder().build_hparams())
        .validate_no_error();
    fix.net.execute_query(fix.conn, "SELECT * FROM empty_table", r).validate_no_error();
}

// any_connection can reconnect after close, even if the stream uses ssl
BOOST_DATA_TEST_CASE(reconnect_after_close_any, any_samples_grid, p0, p1)
{
    netfn_fixture_any fix(p0);
    ssl_mode mode = p1;

    // Connect and use the connection
    results r;
    fix.net.connect(fix.conn, connect_params_builder().ssl(mode).build()).validate_no_error();
    fix.net.execute_query(fix.conn, "SELECT * FROM empty_table", r).validate_no_error();

    // Close
    fix.net.close(fix.conn).validate_no_error();

    // Reopen and use the connection normally
    fix.net.connect(fix.conn, connect_params_builder().ssl(mode).build()).validate_no_error();
    fix.net.execute_query(fix.conn, "SELECT * FROM empty_table", r).validate_no_error();
}

// Old connection can reconnect after handshake failure if the stream is not SSL
BOOST_DATA_TEST_CASE(reconnect_after_handshake_error_connection, connection_samples)
{
    netfn_fixture_connection fix(sample);

    // Error during server handshake
    fix.net.connect(fix.conn, get_tcp_endpoint(), connect_params_builder().database("bad_db").build_hparams())
        .validate_error(
            common_server_errc::er_dbaccess_denied_error,
            "Access denied for user 'integ_user'@'%' to database 'bad_db'"
        );

    // Reopen with correct parameters and use the connection normally
    results r;
    fix.net.connect(fix.conn, get_tcp_endpoint(), connect_params_builder().build_hparams())
        .validate_no_error();
    fix.net.execute_query(fix.conn, "SELECT * FROM empty_table", r).validate_no_error();
}

// any_connection can reconnect after a handshake failure, even if SSL is used
BOOST_DATA_TEST_CASE(reconnect_after_handshake_error_any, any_samples_grid, p0, p1)
{
    netfn_fixture_any fix(p0);
    ssl_mode mode = p1;

    // Error during server handshake
    fix.net.connect(fix.conn, connect_params_builder().ssl(mode).database("bad_db").build())
        .validate_error(
            common_server_errc::er_dbaccess_denied_error,
            "Access denied for user 'integ_user'@'%' to database 'bad_db'"
        );

    // Reopen with correct parameters and use the connection normally
    results r;
    fix.net.connect(fix.conn, connect_params_builder().ssl(mode).build()).validate_no_error();
    fix.net.execute_query(fix.conn, "SELECT * FROM empty_table", r).validate_no_error();
}

// any_connection can reconnect while it's connected
BOOST_DATA_TEST_CASE(reconnect_while_connected, any_samples_grid, p0, p1)
{
    netfn_fixture_any fix(p0);
    ssl_mode mode = p1;

    // Connect and use the connection
    results r;
    fix.net.connect(fix.conn, connect_params_builder().ssl(mode).build()).validate_no_error();
    fix.net.execute_query(fix.conn, "SELECT * FROM empty_table", r).validate_no_error();

    // We can safely connect again
    fix.net.connect(fix.conn, connect_params_builder().ssl(mode).credentials("root", "").build())
        .validate_no_error();

    // We've logged in as root. May be root@% or root@localhost
    fix.net.execute_query(fix.conn, "SELECT SUBSTRING_INDEX(CURRENT_USER(), '@', 1)", r).validate_no_error();
    BOOST_TEST(r.rows() == makerows(1, "root"), per_element());

    // Close succeeds
    fix.net.close(fix.conn).validate_no_error();
}

BOOST_FIXTURE_TEST_CASE(reconnect_after_cancel, any_connection_fixture)
{
    // Setup
    results r;
    connect();

    // Kick an operation that ends up cancelled
    asio::cancellation_signal sig;
    auto netres = conn.async_execute(
        "DO SLEEP(2)",
        r,
        asio::bind_cancellation_slot(sig.slot(), as_netresult)
    );

    // Return to the event loop and emit the signal
    asio::post(asio::bind_executor(ctx, [&]() {
        // Emit the signal
        sig.emit(asio::cancellation_type::terminal);
    }));

    // Wait for the operation to finish
    std::move(netres).validate_error(asio::error::operation_aborted);

    // We can connect again and use the connection
    connect();
    conn.async_execute("SELECT 42", r, as_netresult).validate_no_error();
}

// any_connection can change the stream type used by successive connect calls
BOOST_AUTO_TEST_CASE(change_stream_type)
{
    // Helpers
    const auto tcp_params = connect_params_builder().ssl(ssl_mode::disable).build();
    const auto tcp_ssl_params = connect_params_builder().ssl(ssl_mode::require).build();
    const auto unix_params = connect_params_builder().set_unix().build();

    // Samples
    struct test_case_t
    {
        string_view name;
        const network_functions_any& fns;
        connect_params first_params;
        connect_params second_params;
    };

    std::vector<test_case_t> test_cases{
        {"sync_tcp_tcpssl",  any_samples[0], tcp_params,     tcp_ssl_params},
        {"async_tcp_tcpssl", any_samples[1], tcp_params,     tcp_ssl_params},
        {"async_tcpssl_tcp", any_samples[1], tcp_ssl_params, tcp_params    },
    };
#ifdef BOOST_ASIO_HAS_LOCAL_SOCKETS
    if (get_server_features().unix_sockets)
    {
        test_cases.push_back({"sync_unix_tcpssl", any_samples[0], unix_params, tcp_ssl_params});
        test_cases.push_back({"async_unix_tcpssl", any_samples[1], unix_params, tcp_ssl_params});
        test_cases.push_back({"async_tcpssl_unix", any_samples[1], tcp_ssl_params, unix_params});
        test_cases.push_back({"async_tcp_unix", any_samples[1], tcp_params, unix_params});
    }
#endif

    // Actual test
    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // Setup
            netfn_fixture_any fix(tc.fns);

            // Connect with the first stream type
            fix.net.connect(fix.conn, tc.first_params).validate_no_error();
            fix.net.ping(fix.conn).validate_no_error();

            // Connect with the second stream type
            fix.net.connect(fix.conn, tc.second_params).validate_no_error();
            fix.net.ping(fix.conn).validate_no_error();
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()  // test_reconnect

}  // namespace
