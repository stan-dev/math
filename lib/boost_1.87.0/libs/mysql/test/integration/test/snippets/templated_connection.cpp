//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/defaults.hpp>
#include <boost/mysql/handshake_params.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/tcp.hpp>
#include <boost/mysql/tcp_ssl.hpp>
#include <boost/mysql/unix.hpp>
#include <boost/mysql/with_diagnostics.hpp>

#include <boost/asio/awaitable.hpp>
#include <boost/asio/connect.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/test/unit_test.hpp>

#include "test_common/ci_server.hpp"
#include "test_common/io_context_fixture.hpp"
#include "test_integration/run_coro.hpp"
#include "test_integration/server_features.hpp"
#include "test_integration/snippets/credentials.hpp"

namespace mysql = boost::mysql;
namespace asio = boost::asio;
using namespace mysql::test;

namespace {

BOOST_AUTO_TEST_SUITE(section_templated_connection)

BOOST_AUTO_TEST_CASE(creation)
{
    auto server_hostname = get_hostname();

    //[templated_connection_creation
    // The execution context, required for all I/O operations
    asio::io_context ctx;

    // The SSL context, required for connections that use TLS.
    asio::ssl::context ssl_ctx(asio::ssl::context::tlsv12_client);

    // Construct the connection. The arguments are forwarded
    // to the stream type (asio::ssl::stream<asio::ip::tcp::socket>).
    mysql::tcp_ssl_connection conn(ctx, ssl_ctx);
    //]

    //[templated_connection_connect
    // Resolve the hostname to get a collection of endpoints.
    // default_port_string is MySQL's default port, 3306
    // Hostname resolution may yield more than one host
    asio::ip::tcp::resolver resolver(ctx);
    auto endpoints = resolver.resolve(server_hostname, mysql::default_port_string);

    // Parameters specifying how to perform the MySQL handshake operation.
    // Similar to connect_params, but doesn't contain the server address and is non-owning
    mysql::handshake_params params(
        mysql_username,
        mysql_password,
        "boost_mysql_examples"  // database to use
    );

    // Connect to the server using the first endpoint returned by the resolver
    conn.connect(*endpoints.begin(), params);
    //]

    //[templated_connection_use
    // Issue a query, as you would with any_connection
    mysql::results result;
    conn.execute("SELECT 1", result);
    //]

    //[templated_connection_close
    conn.close();
    //]
}

#ifdef BOOST_ASIO_HAS_LOCAL_SOCKETS
BOOST_TEST_DECORATOR(*run_if(&server_features::unix_sockets))
BOOST_AUTO_TEST_CASE(unix_sockets)
{
    //[templated_connection_unix
    // The execution context, required for all I/O operations
    asio::io_context ctx;

    // A UNIX connection requires only an execution context
    mysql::unix_connection conn(ctx);

    // The socket path where the server is listening
    asio::local::stream_protocol::endpoint ep("/var/run/mysqld/mysqld.sock");

    // MySQL handshake parameters, as in the TCP case.
    mysql::handshake_params params(
        mysql_username,
        mysql_password,
        "boost_mysql_examples"  // database to use
    );

    // Connect to the server
    conn.connect(ep, params);

    // Use the connection normally
    //]
}
#endif

BOOST_AUTO_TEST_CASE(handshake_quit)
{
    auto server_hostname = get_hostname();

    //[templated_connection_handshake_quit
    // The execution context, required for all I/O operations
    asio::io_context ctx;

    // The SSL context, required for connections that use TLS.
    asio::ssl::context ssl_ctx(asio::ssl::context::tlsv12_client);

    // We're using TLS over TCP
    mysql::tcp_ssl_connection conn(ctx, ssl_ctx);

    // Resolve the server hostname into endpoints
    asio::ip::tcp::resolver resolver(ctx);
    auto endpoints = resolver.resolve(server_hostname, mysql::default_port_string);

    // Connect the underlying stream manually.
    // asio::connect tries every endpoint in the passed sequence
    // until one succeeds. any_connection uses this internally.
    // lowest_layer obtains the underlying socket from the ssl::stream
    asio::connect(conn.stream().lowest_layer(), endpoints);

    // Perform MySQL session establishment.
    // This will also perform the TLS handshake, if required.
    mysql::handshake_params params(
        mysql_username,
        mysql_password,
        "boost_mysql_examples"  // database to use
    );
    conn.handshake(params);

    // Use the connection normally
    mysql::results result;
    conn.execute("SELECT 1", result);

    // Terminate the connection. This also performs the TLS shutdown.
    conn.quit();

    // Close the underlying stream.
    // The connection's destructor also closes the socket,
    // but doing it explicitly will throw in case of error.
    conn.stream().lowest_layer().close();
    //]
}

#ifdef BOOST_ASIO_HAS_CO_AWAIT
BOOST_FIXTURE_TEST_CASE(with_diagnostics_, io_context_fixture)
{
    run_coro(ctx, [&]() -> asio::awaitable<void> {
        // Setup
        mysql::tcp_connection conn(ctx);
        asio::ip::tcp::resolver resolv(ctx);
        auto endpoints = co_await resolv.async_resolve(get_hostname(), mysql::default_port_string);
        mysql::handshake_params hparams(mysql_username, mysql_password);
        co_await conn.async_connect(*endpoints.begin(), hparams, mysql::with_diagnostics(asio::deferred));
        mysql::results result;

        //[templated_connection_with_diagnostics
        co_await conn.async_execute("SELECT 1", result, mysql::with_diagnostics(asio::deferred));
        //]
    });
}
#endif

BOOST_AUTO_TEST_SUITE_END()

}  // namespace
