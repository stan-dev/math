//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/connect_params.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/metadata_mode.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/static_execution_state.hpp>

#include <boost/asio/awaitable.hpp>
#include <boost/asio/cancel_after.hpp>
#include <boost/asio/deferred.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/core/span.hpp>
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <exception>
#include <stdexcept>
#include <tuple>

#include "test_common/printing.hpp"
#include "test_unit/test_any_connection.hpp"

using namespace boost::mysql;
namespace asio = boost::asio;
using asio::deferred;

BOOST_AUTO_TEST_SUITE(test_any_connection)

BOOST_AUTO_TEST_CASE(init_ctor)
{
    asio::io_context ctx;
    any_connection c{ctx.get_executor()};
    BOOST_TEST((c.get_executor() == ctx.get_executor()));
}

BOOST_AUTO_TEST_CASE(init_ctor_execution_context)
{
    asio::io_context ctx;
    any_connection c{ctx};
    BOOST_TEST((c.get_executor() == ctx.get_executor()));
}

BOOST_AUTO_TEST_CASE(init_ctor_with_buffer_size)
{
    asio::io_context ctx;
    any_connection_params params;
    params.initial_buffer_size = 512u;
    any_connection c{ctx.get_executor(), params};
    BOOST_TEST((c.get_executor() == ctx.get_executor()));
}

BOOST_AUTO_TEST_CASE(init_ctor_max_buffer_size_eq_size)
{
    asio::io_context ctx;
    any_connection_params params;
    params.initial_buffer_size = 512u;
    params.max_buffer_size = 512u;
    any_connection c{ctx.get_executor(), params};
    BOOST_TEST((c.get_executor() == ctx.get_executor()));
}

BOOST_AUTO_TEST_CASE(init_ctor_error_max_buffer_size)
{
    asio::io_context ctx;
    any_connection_params params;
    params.initial_buffer_size = 513u;
    params.max_buffer_size = 512u;
    BOOST_CHECK_THROW(any_connection(ctx, params), std::invalid_argument);
}

// move ctor
BOOST_AUTO_TEST_CASE(move_ctor)
{
    asio::io_context ctx;
    any_connection c1{ctx.get_executor()};
    any_connection c2{std::move(c1)};
    BOOST_TEST((c2.get_executor() == ctx.get_executor()));
}

// move assign
BOOST_AUTO_TEST_CASE(move_assign_to_moved_from)
{
    asio::io_context ctx;
    any_connection moved_from{ctx.get_executor()};
    any_connection other{std::move(moved_from)};
    any_connection conn{ctx.get_executor()};
    moved_from = std::move(conn);
    BOOST_TEST((moved_from.get_executor() == ctx.get_executor()));
}

BOOST_AUTO_TEST_CASE(move_assign_to_valid)
{
    asio::io_context ctx;
    any_connection c1{ctx.get_executor()};
    any_connection c2{ctx.get_executor()};
    c1 = std::move(c2);
    BOOST_TEST((c1.get_executor() == ctx.get_executor()));
}

BOOST_AUTO_TEST_CASE(set_meta_mode)
{
    asio::io_context ctx;

    // Default metadata mode
    any_connection conn{ctx.get_executor()};
    BOOST_TEST(conn.meta_mode() == metadata_mode::minimal);

    // Setting it causes effect
    conn.set_meta_mode(metadata_mode::full);
    BOOST_TEST(conn.meta_mode() == metadata_mode::full);
}

// spotcheck: deferred compiles even in C++11
void deferred_spotcheck()
{
    asio::io_context ctx;
    auto conn = test::create_test_any_connection(ctx);
    connect_params params;
    diagnostics diag;
    results result;
    execution_state st;
    statement stmt;

    (void)conn.async_connect(params, deferred);
    (void)conn.async_connect(params, diag, deferred);

    (void)conn.async_execute("SELECT 1", result, deferred);
    (void)conn.async_execute("SELECT 1", result, diag, deferred);

    (void)conn.async_start_execution("SELECT 1", st, deferred);
    (void)conn.async_start_execution("SELECT 1", st, diag, deferred);

    (void)conn.async_read_some_rows(st, deferred);
    (void)conn.async_read_some_rows(st, diag, deferred);

#ifdef BOOST_MYSQL_CXX14
    static_execution_state<std::tuple<>> st2;
    (void)conn.async_read_some_rows(st2, boost::span<std::tuple<>>{}, deferred);
    (void)conn.async_read_some_rows(st2, boost::span<std::tuple<>>{}, diag, deferred);
#endif

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

    (void)conn.async_close(deferred);
    (void)conn.async_close(diag, deferred);
}

// Spotcheck: all any_connection functions support default completion tokens
#ifdef BOOST_ASIO_HAS_CO_AWAIT
asio::awaitable<void> spotcheck_default_tokens()
{
    asio::io_context ctx;
    auto conn = test::create_test_any_connection(ctx);
    connect_params params;
    diagnostics diag;
    results result;
    execution_state st;
    statement stmt;
    static_execution_state<std::tuple<>> st2;

    co_await conn.async_connect(params);
    co_await conn.async_connect(params, diag);

    co_await conn.async_execute("SELECT 1", result);
    co_await conn.async_execute("SELECT 1", result, diag);

    co_await conn.async_start_execution("SELECT 1", st);
    co_await conn.async_start_execution("SELECT 1", st, diag);

    co_await conn.async_read_some_rows(st);
    co_await conn.async_read_some_rows(st, diag);

    co_await conn.async_read_some_rows(st2, boost::span<std::tuple<>>{});
    co_await conn.async_read_some_rows(st2, boost::span<std::tuple<>>{}, diag);

    co_await conn.async_read_resultset_head(st);
    co_await conn.async_read_resultset_head(st, diag);

    co_await conn.async_prepare_statement("SELECT 1");
    co_await conn.async_prepare_statement("SELECT 1", diag);

    co_await conn.async_close_statement(stmt);
    co_await conn.async_close_statement(stmt, diag);

    co_await conn.async_reset_connection();
    co_await conn.async_reset_connection(diag);

    co_await conn.async_ping();
    co_await conn.async_ping(diag);

    co_await conn.async_close();
    co_await conn.async_close(diag);
}
#endif

// Spotcheck: any_connection ops support partial tokens,
// and they get passed the correct default token
template <class T, class... SigArgs, class... Rest>
void check_op(asio::deferred_async_operation<void(T, SigArgs...), Rest...>)
{
    static_assert(std::is_same<T, std::exception_ptr>::value, "");
}

void spotcheck_partial_tokens()
{
    asio::io_context ctx;
    auto conn = test::create_test_any_connection(ctx);
    connect_params params;
    diagnostics diag;
    results result;
    execution_state st;
    statement stmt;
    auto tok = asio::cancel_after(std::chrono::seconds(10));

    check_op(conn.async_connect(params, tok));
    check_op(conn.async_connect(params, diag, tok));

    check_op(conn.async_execute("SELECT 1", result, tok));
    check_op(conn.async_execute("SELECT 1", result, diag, tok));

    check_op(conn.async_start_execution("SELECT 1", st, tok));
    check_op(conn.async_start_execution("SELECT 1", st, diag, tok));

    check_op(conn.async_read_some_rows(st, tok));
    check_op(conn.async_read_some_rows(st, diag, tok));

#ifdef BOOST_MYSQL_CXX14
    static_execution_state<std::tuple<>> st2;
    check_op(conn.async_read_some_rows(st2, boost::span<std::tuple<>>{}, tok));
    check_op(conn.async_read_some_rows(st2, boost::span<std::tuple<>>{}, diag, tok));
#endif

    check_op(conn.async_read_resultset_head(st, tok));
    check_op(conn.async_read_resultset_head(st, diag, tok));

    check_op(conn.async_prepare_statement("SELECT 1", tok));
    check_op(conn.async_prepare_statement("SELECT 1", diag, tok));

    check_op(conn.async_close_statement(stmt, tok));
    check_op(conn.async_close_statement(stmt, diag, tok));

    check_op(conn.async_reset_connection(tok));
    check_op(conn.async_reset_connection(diag, tok));

    check_op(conn.async_ping(tok));
    check_op(conn.async_ping(diag, tok));

    check_op(conn.async_close(tok));
    check_op(conn.async_close(diag, tok));
}

BOOST_AUTO_TEST_SUITE_END()
