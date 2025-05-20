//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/buffer_params.hpp>
#include <boost/mysql/connection.hpp>
#include <boost/mysql/handshake_params.hpp>
#include <boost/mysql/metadata_mode.hpp>
#include <boost/mysql/tcp.hpp>
#include <boost/mysql/tcp_ssl.hpp>

#include <boost/asio/deferred.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/strand.hpp>
#include <boost/test/unit_test.hpp>

#include "test_common/printing.hpp"

using namespace boost::mysql;
namespace asio = boost::asio;
using asio::deferred;

BOOST_AUTO_TEST_SUITE(test_connection_)

BOOST_AUTO_TEST_CASE(init_ctor)
{
    asio::io_context ctx;
    tcp_connection c{ctx.get_executor()};
    BOOST_TEST((c.get_executor() == ctx.get_executor()));
    BOOST_TEST((c.stream().get_executor() == ctx.get_executor()));
}

BOOST_AUTO_TEST_CASE(init_ctor_with_buffer_params)
{
    asio::io_context ctx;
    tcp_connection c{buffer_params(4096), ctx.get_executor()};
    BOOST_TEST((c.get_executor() == ctx.get_executor()));
    BOOST_TEST((c.stream().get_executor() == ctx.get_executor()));
}

BOOST_AUTO_TEST_CASE(init_ctor_with_buffer_params_rvalue)
{
    asio::io_context ctx;
    tcp_connection::stream_type sock{ctx.get_executor()};
    tcp_connection c{buffer_params(4096), std::move(sock)};
    BOOST_CHECK_NO_THROW(c.stream().get_executor());
}

// move ctor
BOOST_AUTO_TEST_CASE(move_ctor)
{
    asio::io_context ctx;
    tcp_connection c1{ctx.get_executor()};
    tcp_connection c2{std::move(c1)};
    BOOST_TEST((c2.get_executor() == ctx.get_executor()));
}

// move assign
BOOST_AUTO_TEST_CASE(move_assign_to_moved_from)
{
    asio::io_context ctx;
    tcp_connection moved_from{ctx.get_executor()};
    tcp_connection other{std::move(moved_from)};
    tcp_connection conn{ctx.get_executor()};
    moved_from = std::move(conn);
    BOOST_TEST((moved_from.get_executor() == ctx.get_executor()));
}

BOOST_AUTO_TEST_CASE(move_assign_to_valid)
{
    asio::io_context ctx;
    tcp_connection c1{ctx.get_executor()};
    tcp_connection c2{ctx.get_executor()};
    c1 = std::move(c2);
    BOOST_TEST((c1.get_executor() == ctx.get_executor()));
}

BOOST_AUTO_TEST_CASE(set_meta_mode)
{
    asio::io_context ctx;

    // Default metadata mode
    tcp_connection conn{ctx.get_executor()};
    BOOST_TEST(conn.meta_mode() == metadata_mode::minimal);

    // Setting it causes effect
    conn.set_meta_mode(metadata_mode::full);
    BOOST_TEST(conn.meta_mode() == metadata_mode::full);
}

// rebind_executor
using other_exec = asio::strand<asio::any_io_executor>;
static_assert(
    std::is_same<
        typename tcp_connection::rebind_executor<other_exec>::other,
        connection<asio::basic_stream_socket<asio::ip::tcp, other_exec>>>::value,
    ""
);
static_assert(
    std::is_same<
        typename tcp_ssl_connection::rebind_executor<other_exec>::other,
        connection<asio::ssl::stream<asio::basic_stream_socket<asio::ip::tcp, other_exec>>>>::value,
    ""
);

// Verify that connection doesn't cause compilation errors for a minimal stream type
struct stream_archetype
{
    // Executor
    using executor_type = asio::any_io_executor;
    executor_type get_executor() { return {}; }

    // Reading
    std::size_t read_some(asio::mutable_buffer, error_code&) { return 0; }

    template <class CompletionToken>
    void async_read_some(asio::mutable_buffer, CompletionToken&&)
    {
    }

    // Writing
    std::size_t write_some(asio::const_buffer, error_code&) { return 0; }

    template <class CompletionToken>
    void async_write_some(asio::const_buffer, CompletionToken&&)
    {
    }
};

connection<stream_archetype> archetype_conn;

// spotcheck: deferred compiles even in C++11
void deferred_spotcheck()
{
    asio::io_context ctx;
    tcp_connection conn(ctx);
    handshake_params params("myuser", "mypasswd");
    asio::ip::tcp::endpoint ep;
    diagnostics diag;
    results result;
    execution_state st;
    statement stmt;
    std::string str;
    const std::string const_str;

    (void)conn.async_handshake(params, deferred);
    (void)conn.async_handshake(params, diag, deferred);

    (void)conn.async_execute("SELECT 1", result, deferred);
    (void)conn.async_execute(str, result, deferred);
    (void)conn.async_execute(const_str, result, deferred);
    (void)conn.async_execute("SELECT 1", result, diag, deferred);
    (void)conn.async_execute(str, result, diag, deferred);
    (void)conn.async_execute(const_str, result, diag, deferred);

    (void)conn.async_start_execution("SELECT 1", st, deferred);
    (void)conn.async_start_execution(str, st, deferred);
    (void)conn.async_start_execution(const_str, st, deferred);
    (void)conn.async_start_execution("SELECT 1", st, diag, deferred);
    (void)conn.async_start_execution(str, st, diag, deferred);
    (void)conn.async_start_execution(const_str, st, diag, deferred);

    (void)conn.async_read_some_rows(st, deferred);
    (void)conn.async_read_some_rows(st, diag, deferred);

    (void)conn.async_read_resultset_head(st, deferred);
    (void)conn.async_read_resultset_head(st, diag, deferred);

    (void)conn.async_prepare_statement("SELECT 1", deferred);
    (void)conn.async_prepare_statement("SELECT 1", diag, deferred);

    (void)conn.async_close_statement(stmt, deferred);
    (void)conn.async_close_statement(stmt, diag, deferred);

    (void)conn.async_reset_connection(deferred);
    (void)conn.async_reset_connection(diag, deferred);

    (void)conn.async_ping(deferred);
    (void)conn.async_ping(diag, deferred);

    (void)conn.async_quit(deferred);
    (void)conn.async_quit(diag, deferred);

    (void)conn.async_close(deferred);
    (void)conn.async_close(diag, deferred);
}

BOOST_AUTO_TEST_SUITE_END()  // test_connection
