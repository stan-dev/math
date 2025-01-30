//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/handshake_params.hpp>
#include <boost/mysql/static_execution_state.hpp>
#include <boost/mysql/tcp_ssl.hpp>

#include <boost/asio/awaitable.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/ssl/context.hpp>
#include <boost/asio/use_awaitable.hpp>

#include <tuple>

#ifdef BOOST_ASIO_HAS_CO_AWAIT

namespace boost {
namespace mysql {
namespace test {

// Verify that default completion tokens compile.
// This is just a spotcheck, this function is never called
boost::asio::awaitable<void> test_default_completion_tokens()
{
    using conn_type = boost::asio::use_awaitable_t<>::as_default_on_t<boost::mysql::tcp_ssl_connection>;

    // Connection
    boost::asio::io_context ctx;
    boost::asio::ssl::context ssl_ctx(boost::asio::ssl::context::tls_client);
    conn_type conn{ctx, ssl_ctx};

    // Helpers
    diagnostics diag;
    boost::asio::ip::tcp::endpoint ep{};
    handshake_params params{"user", "pass"};
    results result;
    static_execution_state<std::tuple<>> static_st;
    execution_state st;
    statement stmt;
    std::vector<field_view> stmt_params;

    // Tests
    co_await conn.async_connect(ep, params);
    co_await conn.async_connect(ep, params, diag);

    co_await conn.async_handshake(params);
    co_await conn.async_handshake(params, diag);

    co_await conn.async_prepare_statement("SELECT 1");
    co_await conn.async_prepare_statement("SELECT 1", diag);

    co_await conn.async_execute("SELECT 1", result);
    co_await conn.async_execute("SELECT 1", result, diag);

    co_await conn.async_start_execution("SELECT 1", st);
    co_await conn.async_start_execution("SELECT 1", st, diag);

    co_await conn.async_close_statement(stmt);
    co_await conn.async_close_statement(stmt, diag);

    co_await conn.async_read_resultset_head(st);
    co_await conn.async_read_resultset_head(st, diag);

    co_await conn.async_read_some_rows(st);
    co_await conn.async_read_some_rows(st, diag);

    co_await conn.async_ping();
    co_await conn.async_ping(diag);

    co_await conn.async_quit();
    co_await conn.async_quit(diag);

    co_await conn.async_close();
    co_await conn.async_close(diag);

    co_await conn.async_read_some_rows(static_st, boost::span<std::tuple<>>());
    co_await conn.async_read_some_rows(static_st, boost::span<std::tuple<>>(), diag);
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
