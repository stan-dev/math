//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_SRC_ER_IMPL_COMMON_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_SRC_ER_IMPL_COMMON_HPP

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

#include "test_common/network_result.hpp"
#include "test_integration/er_connection.hpp"
#include "test_integration/er_network_variant.hpp"
#include "test_integration/get_endpoint.hpp"
#include "test_integration/static_rows.hpp"
#include "test_integration/streams.hpp"

namespace boost {
namespace mysql {
namespace test {

// Variants
void add_sync_errc(std::vector<er_network_variant*>&);
void add_sync_exc(std::vector<er_network_variant*>&);
void add_async_callback(std::vector<er_network_variant*>&);
void add_async_coroutines(std::vector<er_network_variant*>&);
void add_async_coroutinescpp20(std::vector<er_network_variant*>&);

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

// Base for implementing er_connection
template <class Stream>
class connection_base : public er_connection
{
public:
    using stream_type = Stream;
    using connection_type = connection<Stream>;
    using stmt_tuple = bound_statement_tuple<std::tuple<field_view, field_view>>;
    using stmt_it = bound_statement_iterator_range<fv_list_it>;

    connection_base(
        boost::asio::any_io_executor executor,
        boost::asio::ssl::context& ssl_ctx,
        er_network_variant& var
    )
        : conn_(create_connection<Stream>(executor, ssl_ctx)), var_(var)
    {
    }

    connection_type& conn() noexcept { return conn_; }

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

private:
    connection_type conn_;
    er_network_variant& var_;
};

// Macros to help implementing er_connection
#ifdef BOOST_MYSQL_CXX14
#define BOOST_MYSQL_TEST_IMPLEMENT_GENERIC_CXX14(prefix)                                                   \
    network_result<void> execute(string_view q, er_connection::static_results_t& result) override          \
    {                                                                                                      \
        return fn_impl<void, const string_view&, er_connection::static_results_t&>(                        \
            &conn_type::prefix##execute,                                                                   \
            q,                                                                                             \
            result                                                                                         \
        );                                                                                                 \
    }                                                                                                      \
    network_result<void> start_execution(string_view q, er_connection::static_state_t& st) override        \
    {                                                                                                      \
        return fn_impl<void, const string_view&, er_connection::static_state_t&>(                          \
            &conn_type::prefix##start_execution,                                                           \
            q,                                                                                             \
            st                                                                                             \
        );                                                                                                 \
    }                                                                                                      \
    network_result<void> read_resultset_head(er_connection::static_state_t& st) override                   \
    {                                                                                                      \
        return fn_impl<void, er_connection::static_state_t&>(&conn_type::prefix##read_resultset_head, st); \
    }                                                                                                      \
    network_result<std::size_t> read_some_rows(                                                            \
        er_connection::static_state_t& st,                                                                 \
        boost::span<row_multifield> storage                                                                \
    ) override                                                                                             \
    {                                                                                                      \
        return fn_impl<std::size_t, er_connection::static_state_t&, boost::span<row_multifield>>(          \
            &conn_type::prefix##read_some_rows,                                                            \
            st,                                                                                            \
            storage                                                                                        \
        );                                                                                                 \
    }                                                                                                      \
    network_result<std::size_t> read_some_rows(                                                            \
        er_connection::static_state_t& st,                                                                 \
        boost::span<row_2fields> storage                                                                   \
    ) override                                                                                             \
    {                                                                                                      \
        return fn_impl<std::size_t, er_connection::static_state_t&, boost::span<row_2fields>>(             \
            &conn_type::prefix##read_some_rows,                                                            \
            st,                                                                                            \
            storage                                                                                        \
        );                                                                                                 \
    }
#else
#define BOOST_MYSQL_TEST_IMPLEMENT_GENERIC_CXX14(prefix)
#endif

#define BOOST_MYSQL_TEST_IMPLEMENT_GENERIC(prefix)                                                           \
    using connection_base<Stream>::connection_base;                                                          \
    network_result<void> connect(const handshake_params& params) override                                    \
    {                                                                                                        \
        return fn_impl<                                                                                      \
            void,                                                                                            \
            const typename Stream::lowest_layer_type::endpoint_type&,                                        \
            const handshake_params&>(&conn_type::prefix##connect, get_endpoint<Stream>(), params);           \
    }                                                                                                        \
    network_result<void> handshake(const handshake_params& params) override                                  \
    {                                                                                                        \
        return fn_impl<void, const handshake_params&>(&conn_type::prefix##handshake, params);                \
    }                                                                                                        \
    network_result<void> query(string_view query, results& result) override                                  \
    {                                                                                                        \
        return fn_impl<void, string_view, results&>(&conn_type::prefix##query, query, result);               \
    }                                                                                                        \
    network_result<void> start_query(string_view query, execution_state& st) override                        \
    {                                                                                                        \
        return fn_impl<void, string_view, execution_state&>(&conn_type::prefix##start_query, query, st);     \
    }                                                                                                        \
    network_result<statement> prepare_statement(string_view stmt_sql) override                               \
    {                                                                                                        \
        return fn_impl<statement, string_view>(&conn_type::prefix##prepare_statement, stmt_sql);             \
    }                                                                                                        \
    network_result<void>                                                                                     \
    execute_statement(const statement& stmt, field_view param1, field_view param2, results& result) override \
    {                                                                                                        \
        return fn_impl<void, const statement&, const std::tuple<field_view, field_view>&, results&>(         \
            &conn_type::prefix##execute_statement,                                                           \
            stmt,                                                                                            \
            std::make_tuple(param1, param2),                                                                 \
            result                                                                                           \
        );                                                                                                   \
    }                                                                                                        \
    network_result<void> start_statement_execution(                                                          \
        const statement& stmt,                                                                               \
        field_view param1,                                                                                   \
        field_view param2,                                                                                   \
        execution_state& st                                                                                  \
    ) override                                                                                               \
    {                                                                                                        \
        return fn_impl<void, const statement&, const std::tuple<field_view, field_view>&, execution_state&>( \
            &conn_type::prefix##start_statement_execution,                                                   \
            stmt,                                                                                            \
            std::make_tuple(param1, param2),                                                                 \
            st                                                                                               \
        );                                                                                                   \
    }                                                                                                        \
    network_result<void> start_statement_execution(                                                          \
        const statement& stmt,                                                                               \
        fv_list_it params_first,                                                                             \
        fv_list_it params_last,                                                                              \
        execution_state& st                                                                                  \
    ) override                                                                                               \
    {                                                                                                        \
        return fn_impl<void, const statement&, fv_list_it, fv_list_it, execution_state&>(                    \
            &conn_type::prefix##start_statement_execution,                                                   \
            stmt,                                                                                            \
            params_first,                                                                                    \
            params_last,                                                                                     \
            st                                                                                               \
        );                                                                                                   \
    }                                                                                                        \
    network_result<void> execute(string_view query, results& result) override                                \
    {                                                                                                        \
        return fn_impl<void, const string_view&, results&>(&conn_type::prefix##execute, query, result);      \
    }                                                                                                        \
    network_result<void> execute(                                                                            \
        bound_statement_tuple<std::tuple<field_view, field_view>> req,                                       \
        results& result                                                                                      \
    ) override                                                                                               \
    {                                                                                                        \
        return fn_impl<void, const typename base_type::stmt_tuple&, results&>(                               \
            &conn_type::prefix##execute,                                                                     \
            req,                                                                                             \
            result                                                                                           \
        );                                                                                                   \
    }                                                                                                        \
    network_result<void> execute(bound_statement_iterator_range<fv_list_it> req, results& result) override   \
    {                                                                                                        \
        return fn_impl<void, const typename base_type::stmt_it&, results&>(                                  \
            &conn_type::prefix##execute,                                                                     \
            req,                                                                                             \
            result                                                                                           \
        );                                                                                                   \
    }                                                                                                        \
    network_result<void> start_execution(string_view query, execution_state& st) override                    \
    {                                                                                                        \
        return fn_impl<void, const string_view&, execution_state&>(                                          \
            &conn_type::prefix##start_execution,                                                             \
            query,                                                                                           \
            st                                                                                               \
        );                                                                                                   \
    }                                                                                                        \
    network_result<void> start_execution(                                                                    \
        bound_statement_tuple<std::tuple<field_view, field_view>> req,                                       \
        execution_state& st                                                                                  \
    ) override                                                                                               \
    {                                                                                                        \
        return fn_impl<void, const typename base_type::stmt_tuple&, execution_state&>(                       \
            &conn_type::prefix##start_execution,                                                             \
            req,                                                                                             \
            st                                                                                               \
        );                                                                                                   \
    }                                                                                                        \
    network_result<void> start_execution(                                                                    \
        bound_statement_iterator_range<fv_list_it> req,                                                      \
        execution_state& st                                                                                  \
    ) override                                                                                               \
    {                                                                                                        \
        return fn_impl<void, const typename base_type::stmt_it&, execution_state&>(                          \
            &conn_type::prefix##start_execution,                                                             \
            req,                                                                                             \
            st                                                                                               \
        );                                                                                                   \
    }                                                                                                        \
    network_result<void> close_statement(statement& stmt) override                                           \
    {                                                                                                        \
        return fn_impl<void, const statement&>(&conn_type::prefix##close_statement, stmt);                   \
    }                                                                                                        \
    network_result<void> read_resultset_head(execution_state& st) override                                   \
    {                                                                                                        \
        return fn_impl<void, execution_state&>(&conn_type::prefix##read_resultset_head, st);                 \
    }                                                                                                        \
    network_result<rows_view> read_some_rows(execution_state& st) override                                   \
    {                                                                                                        \
        return fn_impl<rows_view, execution_state&>(&conn_type::prefix##read_some_rows, st);                 \
    }                                                                                                        \
    network_result<void> ping() override { return fn_impl<void>(&conn_type::prefix##ping); }                 \
    network_result<void> reset_connection() override                                                         \
    {                                                                                                        \
        return fn_impl<void>(&conn_type::prefix##reset_connection);                                          \
    }                                                                                                        \
    network_result<void> quit() override { return fn_impl<void>(&conn_type::prefix##quit); }                 \
    network_result<void> close() override { return fn_impl<void>(&conn_type::prefix##close); }               \
    BOOST_MYSQL_TEST_IMPLEMENT_GENERIC_CXX14(prefix)

// Use these
#define BOOST_MYSQL_TEST_IMPLEMENT_SYNC() BOOST_MYSQL_TEST_IMPLEMENT_GENERIC()
#define BOOST_MYSQL_TEST_IMPLEMENT_ASYNC() BOOST_MYSQL_TEST_IMPLEMENT_GENERIC(async_)

// Implementation for er_network_variant
template <class ErConnection>
class er_network_variant_impl : public er_network_variant
{
    using stream_type = typename ErConnection::stream_type;

public:
    bool supports_ssl() const override { return ::boost::mysql::test::supports_ssl<stream_type>(); }
    bool is_unix_socket() const override { return ::boost::mysql::test::is_unix_socket<stream_type>(); }
    const char* stream_name() const override { return ::boost::mysql::test::get_stream_name<stream_type>(); }
    const char* variant_name() const override { return ErConnection::name(); }
    er_connection_ptr create_connection(boost::asio::any_io_executor ex, boost::asio::ssl::context& ssl_ctx)
        override
    {
        return er_connection_ptr(new ErConnection(std::move(ex), ssl_ctx, *this));
    }
};

// Helpers
inline void rethrow_on_failure(std::exception_ptr ptr)
{
    if (ptr)
    {
        std::rethrow_exception(ptr);
    }
}

template <class ErConnection>
void add_variant(std::vector<er_network_variant*>& output)
{
    static er_network_variant_impl<ErConnection> obj;
    output.push_back(&obj);
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
