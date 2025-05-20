//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/asio/ip/tcp.hpp>

#include <ostream>
#include <utility>

#include "test_integration/spotchecks_helpers.hpp"

using namespace boost::mysql::test;

// Network functions
#define BOOST_MYSQL_MAKE_NETFN(conn, netm, fn, i)              \
    (i == 0   ? netmakers::netm::sync_errc(&conn::fn)          \
     : i == 1 ? netmakers::netm::sync_exc(&conn::fn)           \
     : i == 2 ? netmakers::netm::async_diag(&conn::async_##fn) \
              : netmakers::netm::async_nodiag(&conn::async_##fn))

static constexpr const char* fn_names[] = {"sync_errc", "sync_exc", "async_diag", "async_nodiag"};

std::vector<network_functions_connection> boost::mysql::test::network_functions_connection::all()
{
    std::vector<network_functions_connection> res;

    for (std::size_t i = 0; i < 4; ++i)
    {
        // connect_stream doesn't involve diagnostics
        auto fn_connect_stream = (i == 0 || i == 1) ? netmakers::connect_stream::sync_errc_nodiag(
                                                          &asio::ip::tcp::socket::connect
                                                      )
                                                    : netmakers::connect_stream::async_nodiag(
                                                          &asio::ip::tcp::socket::async_connect
                                                      );

        res.push_back({
            fn_names[i],
            BOOST_MYSQL_MAKE_NETFN(tcp_connection, prepare_statement, prepare_statement, i),
            BOOST_MYSQL_MAKE_NETFN(tcp_connection, execute_query, execute, i),
            BOOST_MYSQL_MAKE_NETFN(tcp_connection, execute_statement, execute, i),
            BOOST_MYSQL_MAKE_NETFN(tcp_connection, start_execution, start_execution, i),
            BOOST_MYSQL_MAKE_NETFN(tcp_connection, close_statement, close_statement, i),
            BOOST_MYSQL_MAKE_NETFN(tcp_connection, read_resultset_head, read_resultset_head, i),
            BOOST_MYSQL_MAKE_NETFN(tcp_connection, read_some_rows, read_some_rows, i),
            BOOST_MYSQL_MAKE_NETFN(tcp_connection, ping, ping, i),
            BOOST_MYSQL_MAKE_NETFN(tcp_connection, reset_connection, reset_connection, i),
            BOOST_MYSQL_MAKE_NETFN(tcp_connection, close, close, i),
#ifdef BOOST_MYSQL_CXX14
            BOOST_MYSQL_MAKE_NETFN(tcp_connection, execute_static, execute, i),
            BOOST_MYSQL_MAKE_NETFN(tcp_connection, start_execution_static, start_execution, i),
            BOOST_MYSQL_MAKE_NETFN(tcp_connection, read_resultset_head_static, read_resultset_head, i),
            BOOST_MYSQL_MAKE_NETFN(tcp_connection, read_some_rows_static_1, read_some_rows, i),
            BOOST_MYSQL_MAKE_NETFN(tcp_connection, read_some_rows_static_2, read_some_rows, i),
#endif
            std::move(fn_connect_stream),
            BOOST_MYSQL_MAKE_NETFN(tcp_connection, handshake, handshake, i),
            BOOST_MYSQL_MAKE_NETFN(tcp_connection, connect, connect, i),
            BOOST_MYSQL_MAKE_NETFN(tcp_connection, quit, quit, i),
        });
    }

    return res;
}

std::vector<network_functions_connection> boost::mysql::test::network_functions_connection::sync_and_async()
{
    auto all_fns = all();
    return {std::move(all_fns[0]), std::move(all_fns[2])};
}

std::ostream& boost::mysql::test::operator<<(std::ostream& os, const network_functions_connection& v)
{
    return os << v.name;
}

std::vector<network_functions_any> boost::mysql::test::network_functions_any::all()
{
    std::vector<network_functions_any> res;

    for (std::size_t i = 0; i < 4; ++i)
    {
        res.push_back({
            fn_names[i],
            BOOST_MYSQL_MAKE_NETFN(any_connection, prepare_statement, prepare_statement, i),
            BOOST_MYSQL_MAKE_NETFN(any_connection, execute_query, execute, i),
            BOOST_MYSQL_MAKE_NETFN(any_connection, execute_statement, execute, i),
            BOOST_MYSQL_MAKE_NETFN(any_connection, start_execution, start_execution, i),
            BOOST_MYSQL_MAKE_NETFN(any_connection, close_statement, close_statement, i),
            BOOST_MYSQL_MAKE_NETFN(any_connection, read_resultset_head, read_resultset_head, i),
            BOOST_MYSQL_MAKE_NETFN(any_connection, read_some_rows, read_some_rows, i),
            BOOST_MYSQL_MAKE_NETFN(any_connection, ping, ping, i),
            BOOST_MYSQL_MAKE_NETFN(any_connection, reset_connection, reset_connection, i),
            BOOST_MYSQL_MAKE_NETFN(any_connection, close, close, i),
#ifdef BOOST_MYSQL_CXX14
            BOOST_MYSQL_MAKE_NETFN(any_connection, execute_static, execute, i),
            BOOST_MYSQL_MAKE_NETFN(any_connection, start_execution_static, start_execution, i),
            BOOST_MYSQL_MAKE_NETFN(any_connection, read_resultset_head_static, read_resultset_head, i),
            BOOST_MYSQL_MAKE_NETFN(any_connection, read_some_rows_static_1, read_some_rows, i),
            BOOST_MYSQL_MAKE_NETFN(any_connection, read_some_rows_static_2, read_some_rows, i),
#endif
            BOOST_MYSQL_MAKE_NETFN(any_connection, connect, connect, i),
            BOOST_MYSQL_MAKE_NETFN(any_connection, set_character_set, set_character_set, i),
            BOOST_MYSQL_MAKE_NETFN(any_connection, run_pipeline, run_pipeline, i),
        });
    }

    return res;
}

std::vector<network_functions_any> boost::mysql::test::network_functions_any::sync_and_async()
{
    auto all_fns = all();
    return {std::move(all_fns[0]), std::move(all_fns[2])};
}

std::ostream& boost::mysql::test::operator<<(std::ostream& os, const network_functions_any& v)
{
    return os << v.name;
}