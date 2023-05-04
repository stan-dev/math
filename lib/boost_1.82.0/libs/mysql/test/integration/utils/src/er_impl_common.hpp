//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_UTILS_SRC_ER_IMPL_COMMON_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_UTILS_SRC_ER_IMPL_COMMON_HPP

#include <boost/mysql/connection.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/handshake_params.hpp>
#include <boost/mysql/metadata_mode.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/rows_view.hpp>
#include <boost/mysql/statement.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/asio/any_io_executor.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/ssl/context.hpp>

#include <type_traits>

#include "er_connection.hpp"
#include "er_network_variant.hpp"
#include "get_endpoint.hpp"
#include "streams.hpp"

namespace boost {
namespace mysql {
namespace test {

// Variants
void add_sync_errc(std::vector<er_network_variant*>&);
void add_sync_exc(std::vector<er_network_variant*>&);
void add_async_callback(std::vector<er_network_variant*>&);
void add_async_coroutines(std::vector<er_network_variant*>&);
void add_async_coroutinescpp20(std::vector<er_network_variant*>&);

// Function table
template <class Stream>
struct function_table
{
    using conn_type = connection<Stream>;

    using connect_sig = network_result<
        void>(conn_type&, const typename Stream::lowest_layer_type::endpoint_type&, const handshake_params&);
    using handshake_sig = network_result<void>(conn_type&, const handshake_params&);
    using query_sig = network_result<void>(conn_type&, string_view, results&);
    using start_query_sig = network_result<void>(conn_type&, string_view, execution_state&);
    using prepare_statement_sig = network_result<statement>(conn_type&, string_view);
    using execute_stmt_sig = network_result<
        void>(conn_type&, const statement&, const std::tuple<field_view, field_view>&, results&);
    using start_stmt_execution_tuple_sig = network_result<
        void>(conn_type&, const statement&, const std::tuple<field_view, field_view>&, execution_state&);
    using start_stmt_execution_it_sig =
        network_result<void>(conn_type&, const statement&, fv_list_it, fv_list_it, execution_state&);
    using close_stmt_sig = network_result<void>(conn_type&, const statement&);
    using read_some_rows_sig = network_result<rows_view>(conn_type&, execution_state&);
    using ping_sig = network_result<void>(conn_type&);
    using quit_sig = network_result<void>(conn_type&);
    using close_sig = network_result<void>(conn_type&);

    std::function<connect_sig> connect;
    std::function<handshake_sig> handshake;
    std::function<query_sig> query;
    std::function<start_query_sig> start_query;
    std::function<prepare_statement_sig> prepare_statement;
    std::function<execute_stmt_sig> execute_stmt;
    std::function<start_stmt_execution_tuple_sig> start_stmt_execution_tuple;
    std::function<start_stmt_execution_it_sig> start_stmt_execution_it;
    std::function<close_stmt_sig> close_stmt;
    std::function<read_some_rows_sig> read_some_rows;
    std::function<ping_sig> ping;
    std::function<quit_sig> quit;
    std::function<close_sig> close;
};

// Note: Netmaker should be a struct with a
// public template type called "type", accepting the signature
// arguments as template parameters, and having a static call that
// takes function pointers and returns a type-erased function
// Generate the function table given a sync maker
template <class Stream, class Netmaker>
function_table<Stream> create_sync_table()
{
    using conn_type = connection<Stream>;
    using table_t = function_table<Stream>;

    // clang-format off
    return function_table<Stream>{
        Netmaker::template type<typename table_t::connect_sig>::call(&conn_type::connect),
        Netmaker::template type<typename table_t::handshake_sig>::call(&conn_type::handshake),
        Netmaker::template type<typename table_t::query_sig>::call(&conn_type::query),
        Netmaker::template type<typename table_t::start_query_sig>::call(&conn_type::start_query),
        Netmaker::template type<typename table_t::prepare_statement_sig>::call(&conn_type::prepare_statement),
        Netmaker::template type<typename table_t::execute_stmt_sig>::call(&conn_type::execute_statement),
        Netmaker::template type<typename table_t::start_stmt_execution_tuple_sig>::call(&conn_type::start_statement_execution),
        Netmaker::template type<typename table_t::start_stmt_execution_it_sig>::call(&conn_type::start_statement_execution),
        Netmaker::template type<typename table_t::close_stmt_sig>::call(&conn_type::close_statement),
        Netmaker::template type<typename table_t::read_some_rows_sig>::call(&conn_type::read_some_rows),
        Netmaker::template type<typename table_t::ping_sig>::call(&conn_type::ping),
        Netmaker::template type<typename table_t::quit_sig>::call(&conn_type::quit),
        Netmaker::template type<typename table_t::close_sig>::call(&conn_type::close),
    };
    // clang-format on
}

// Generate the function table given an async maker
template <class Stream, class Netmaker>
function_table<Stream> create_async_table()
{
    using conn_type = connection<Stream>;
    using table_t = function_table<Stream>;

    // clang-format off
    return function_table<Stream>{
        Netmaker::template type<typename table_t::connect_sig>::call(&conn_type::async_connect),
        Netmaker::template type<typename table_t::handshake_sig>::call(&conn_type::async_handshake),
        Netmaker::template type<typename table_t::query_sig>::call(&conn_type::async_query),
        Netmaker::template type<typename table_t::start_query_sig>::call(&conn_type::async_start_query),
        Netmaker::template type<typename table_t::prepare_statement_sig>::call(&conn_type::async_prepare_statement),
        Netmaker::template type<typename table_t::execute_stmt_sig>::call(&conn_type::async_execute_statement),
        Netmaker::template type<typename table_t::start_stmt_execution_tuple_sig>::call(&conn_type::async_start_statement_execution),
        Netmaker::template type<typename table_t::start_stmt_execution_it_sig>::call(&conn_type::async_start_statement_execution),
        Netmaker::template type<typename table_t::close_stmt_sig>::call(&conn_type::async_close_statement),
        Netmaker::template type<typename table_t::read_some_rows_sig>::call(&conn_type::async_read_some_rows),
        Netmaker::template type<typename table_t::ping_sig>::call(&conn_type::async_ping),
        Netmaker::template type<typename table_t::quit_sig>::call(&conn_type::async_quit),
        Netmaker::template type<typename table_t::close_sig>::call(&conn_type::async_close),
    };
    // clang-format on
}

// Helpers
template <class Stream>
connection<Stream> create_connection_impl(
    boost::asio::any_io_executor executor,
    boost::asio::ssl::context& ssl_ctx,
    std::true_type  // is ssl stream
)
{
    return connection<Stream>(executor, ssl_ctx);
}

template <class Stream>
connection<Stream> create_connection_impl(
    boost::asio::any_io_executor executor,
    boost::asio::ssl::context&,
    std::false_type  // is ssl stream
)
{
    return connection<Stream>(executor);
}

template <class Stream>
connection<Stream> create_connection(
    boost::asio::any_io_executor executor,
    boost::asio::ssl::context& ssl_ctx
)
{
    return create_connection_impl<Stream>(
        executor,
        ssl_ctx,
        std::integral_constant<bool, supports_ssl<Stream>()>()
    );
}

// Implementation for er_connection
template <class Stream>
class er_connection_impl : public er_connection
{
    using conn_type = connection<Stream>;

    conn_type conn_;
    er_network_variant& var_;
    const function_table<Stream>& table_;

public:
    er_connection_impl(
        boost::asio::any_io_executor executor,
        boost::asio::ssl::context& ssl_ctx,
        er_network_variant& var,
        const function_table<Stream>& table
    )
        : conn_(create_connection<Stream>(executor, ssl_ctx)), var_(var), table_(table)
    {
    }

    bool uses_ssl() const override { return conn_.uses_ssl(); }
    bool is_open() const override { return conn_.stream().lowest_layer().is_open(); }
    void set_metadata_mode(metadata_mode v) override { conn_.set_meta_mode(v); }
    void physical_connect() override { conn_.stream().lowest_layer().connect(get_endpoint<Stream>()); }
    void sync_close() noexcept override
    {
        try
        {
            conn_.close();
        }
        catch (...)
        {
        }
    }
    er_network_variant& variant() const override { return var_; }

    network_result<void> connect(const handshake_params& params) override
    {
        return table_.connect(conn_, get_endpoint<Stream>(), params);
    }
    network_result<void> handshake(const handshake_params& params) override
    {
        return table_.handshake(conn_, params);
    }
    network_result<void> query(string_view query, results& result) override
    {
        return table_.query(conn_, query, result);
    }
    network_result<void> start_query(string_view query, execution_state& st) override
    {
        return table_.start_query(conn_, query, st);
    }
    network_result<statement> prepare_statement(string_view stmt_sql) override
    {
        return table_.prepare_statement(conn_, stmt_sql);
    }
    network_result<void> execute_statement(
        const statement& stmt,
        field_view param1,
        field_view param2,
        results& result
    ) override
    {
        return table_.execute_stmt(conn_, stmt, std::make_tuple(param1, param2), result);
    }
    network_result<void> start_statement_execution(
        const statement& stmt,
        field_view param1,
        field_view param2,
        execution_state& st
    ) override
    {
        return table_.start_stmt_execution_tuple(conn_, stmt, std::make_tuple(param1, param2), st);
    }
    network_result<void> start_statement_execution(
        const statement& stmt,
        fv_list_it params_first,
        fv_list_it params_last,
        execution_state& st
    ) override
    {
        return table_.start_stmt_execution_it(conn_, stmt, params_first, params_last, st);
    }
    network_result<void> close_statement(statement& stmt) override { return table_.close_stmt(conn_, stmt); }
    network_result<rows_view> read_some_rows(execution_state& st) override
    {
        return table_.read_some_rows(conn_, st);
    }
    network_result<void> ping() override { return table_.ping(conn_); }
    network_result<void> quit() override { return table_.quit(conn_); }
    network_result<void> close() override { return table_.close(conn_); }
};

template <class Stream>
class er_network_variant_impl : public er_network_variant
{
    function_table<Stream> table_;
    const char* variant_name_;

    using conn_type = er_connection_impl<Stream>;

public:
    er_network_variant_impl(function_table<Stream>&& table, const char* name)
        : table_(std::move(table)), variant_name_(name)
    {
    }

    bool supports_ssl() const override { return ::boost::mysql::test::supports_ssl<Stream>(); }
    bool is_unix_socket() const override { return ::boost::mysql::test::is_unix_socket<Stream>(); }
    const char* stream_name() const override { return ::boost::mysql::test::get_stream_name<Stream>(); }
    const char* variant_name() const override { return variant_name_; }
    er_connection_ptr create_connection(boost::asio::any_io_executor ex, boost::asio::ssl::context& ssl_ctx)
        override
    {
        return er_connection_ptr(new conn_type(ex, ssl_ctx, *this, table_));
    }
};

template <class Stream, class Netmaker>
er_network_variant_impl<Stream> create_sync_variant()
{
    return er_network_variant_impl<Stream>(create_sync_table<Stream, Netmaker>(), Netmaker::name());
}

template <class Stream, class Netmaker>
er_network_variant_impl<Stream> create_async_variant()
{
    return er_network_variant_impl<Stream>(create_async_table<Stream, Netmaker>(), Netmaker::name());
}

// Helpers
inline boost::asio::io_context& get_context(boost::asio::any_io_executor ex) noexcept
{
    return static_cast<boost::asio::io_context&>(ex.context());
}

inline void run_until_completion(boost::asio::any_io_executor ex)
{
    auto& ctx = get_context(ex);
    ctx.restart();
    ctx.run();
}

inline void rethrow_on_failure(std::exception_ptr ptr)
{
    if (ptr)
    {
        std::rethrow_exception(ptr);
    }
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
