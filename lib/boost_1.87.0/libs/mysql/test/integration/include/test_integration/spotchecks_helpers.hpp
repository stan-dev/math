//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_SPOTCHECKS_HELPERS_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_SPOTCHECKS_HELPERS_HPP

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/character_set.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/pipeline.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/statement.hpp>
#include <boost/mysql/static_execution_state.hpp>
#include <boost/mysql/static_results.hpp>
#include <boost/mysql/tcp.hpp>

#include <boost/asio/ip/tcp.hpp>
#include <boost/core/span.hpp>

#include <iosfwd>

#include "test_common/io_context_fixture.hpp"
#include "test_common/netfun_maker.hpp"
#include "test_integration/connect_params_builder.hpp"
#include "test_integration/static_rows.hpp"
#include "test_integration/tcp_connection_fixture.hpp"

namespace boost {
namespace mysql {
namespace test {

#ifdef BOOST_MYSQL_CXX14
using static_results_t = static_results<row_multifield, row_2fields, empty>;
using static_state_t = static_execution_state<row_multifield, row_2fields, empty>;
#endif

// Netmakers
template <class Conn>
struct netmakers_common
{
    using prepare_statement = netfun_maker<statement, Conn, string_view>;
    using execute_query = netfun_maker<void, Conn, const string_view&, results&>;
    using execute_statement = netfun_maker<
        void,
        Conn,
        const bound_statement_tuple<std::tuple<int, int>>&,
        results&>;
    using start_execution = netfun_maker<void, Conn, const string_view&, execution_state&>;
    using close_statement = netfun_maker<void, Conn, const statement&>;
    using read_resultset_head = netfun_maker<void, Conn, execution_state&>;
    using read_some_rows = netfun_maker<rows_view, Conn, execution_state&>;
    using ping = netfun_maker<void, Conn>;
    using reset_connection = netfun_maker<void, Conn>;
    using close = netfun_maker<void, Conn>;

#ifdef BOOST_MYSQL_CXX14
    using execute_static = netfun_maker<void, Conn, const string_view&, static_results_t&>;
    using start_execution_static = netfun_maker<void, Conn, const string_view&, static_state_t&>;
    using read_resultset_head_static = netfun_maker<void, Conn, static_state_t&>;
    using read_some_rows_static_1 = netfun_maker<
        std::size_t,
        Conn,
        static_state_t&,
        boost::span<row_multifield>>;
    using read_some_rows_static_2 = netfun_maker<
        std::size_t,
        Conn,
        static_state_t&,
        boost::span<row_2fields>>;
#endif
};

struct netmakers_connection : netmakers_common<tcp_connection>
{
    using connect_stream = netfun_maker<void, asio::ip::tcp::socket, const asio::ip::tcp::endpoint&>;
    using handshake = netfun_maker<void, tcp_connection, const handshake_params&>;
    using connect = netfun_maker<
        void,
        tcp_connection,
        const asio::ip::tcp::endpoint&,
        const handshake_params&>;
    using quit = netfun_maker<void, tcp_connection>;
};

struct netmakers_any : netmakers_common<any_connection>
{
    using connect = netfun_maker<void, any_connection, const connect_params&>;
    using set_character_set = netfun_maker<void, any_connection, const character_set&>;
    using run_pipeline = netfun_maker<
        void,
        any_connection,
        const pipeline_request&,
        std::vector<stage_response>&>;
};

struct network_functions_connection
{
    using netmakers = netmakers_connection;

    string_view name;
    netmakers::prepare_statement::signature prepare_statement;
    netmakers::execute_query::signature execute_query;
    netmakers::execute_statement::signature execute_statement;
    netmakers::start_execution::signature start_execution;
    netmakers::close_statement::signature close_statement;
    netmakers::read_resultset_head::signature read_resultset_head;
    netmakers::read_some_rows::signature read_some_rows;
    netmakers::ping::signature ping;
    netmakers::reset_connection::signature reset_connection;
    netmakers::close::signature close;
#ifdef BOOST_MYSQL_CXX14
    netmakers::execute_static::signature execute_static;
    netmakers::start_execution_static::signature start_execution_static;
    netmakers::read_resultset_head_static::signature read_resultset_head_static;
    netmakers::read_some_rows_static_1::signature read_some_rows_static_1;
    netmakers::read_some_rows_static_2::signature read_some_rows_static_2;
#endif
    netmakers::connect_stream::signature connect_stream;
    netmakers::handshake::signature handshake;
    netmakers::connect::signature connect;
    netmakers::quit::signature quit;

    static std::vector<network_functions_connection> all();
    static std::vector<network_functions_connection> sync_and_async();
};
std::ostream& operator<<(std::ostream& os, const network_functions_connection& v);

struct network_functions_any
{
    using netmakers = netmakers_any;

    string_view name;
    netmakers::prepare_statement::signature prepare_statement;
    netmakers::execute_query::signature execute_query;
    netmakers::execute_statement::signature execute_statement;
    netmakers::start_execution::signature start_execution;
    netmakers::close_statement::signature close_statement;
    netmakers::read_resultset_head::signature read_resultset_head;
    netmakers::read_some_rows::signature read_some_rows;
    netmakers::ping::signature ping;
    netmakers::reset_connection::signature reset_connection;
    netmakers::close::signature close;
#ifdef BOOST_MYSQL_CXX14
    netmakers::execute_static::signature execute_static;
    netmakers::start_execution_static::signature start_execution_static;
    netmakers::read_resultset_head_static::signature read_resultset_head_static;
    netmakers::read_some_rows_static_1::signature read_some_rows_static_1;
    netmakers::read_some_rows_static_2::signature read_some_rows_static_2;
#endif
    netmakers::connect::signature connect;
    netmakers::set_character_set::signature set_character_set;
    netmakers::run_pipeline::signature run_pipeline;

    static std::vector<network_functions_any> all();
    static std::vector<network_functions_any> sync_and_async();
};
std::ostream& operator<<(std::ostream& os, const network_functions_any& v);

// Fixtures. Like tcp_connection_fixture and any_connection_fixture,
// but using a network_functions value
struct netfn_fixture_connection : io_context_fixture
{
    tcp_connection conn{ctx};
    network_functions_connection net;

    netfn_fixture_connection(network_functions_connection n) : net(std::move(n))
    {
        conn.set_meta_mode(metadata_mode::full);
    }
    ~netfn_fixture_connection() { net.close(conn).validate_no_error(); }

    void connect(connect_params_builder conn_params = {})
    {
        net.connect(conn, get_tcp_endpoint(), conn_params.build_hparams()).validate_no_error();
    }
};

struct netfn_fixture_any : io_context_fixture
{
    any_connection conn{ctx};
    network_functions_any net;

    netfn_fixture_any(network_functions_any n) : net(std::move(n))
    {
        conn.set_meta_mode(metadata_mode::full);
    }
    ~netfn_fixture_any() { net.close(conn).validate_no_error(); }

    void connect(connect_params_builder conn_params = {})
    {
        net.connect(conn, conn_params.build()).validate_no_error();
    }
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
