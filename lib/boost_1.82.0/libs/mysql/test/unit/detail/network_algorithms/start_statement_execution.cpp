//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/blob.hpp>
#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/date.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/statement.hpp>

#include <boost/mysql/detail/auxiliar/access_fwd.hpp>
#include <boost/mysql/detail/network_algorithms/start_statement_execution.hpp>
#include <boost/mysql/detail/protocol/constants.hpp>
#include <boost/mysql/detail/protocol/prepared_statement_messages.hpp>
#include <boost/mysql/detail/protocol/resultset_encoding.hpp>

#include <boost/asio/io_context.hpp>
#include <boost/asio/use_awaitable.hpp>
#include <boost/test/unit_test.hpp>

#include <iterator>
#include <tuple>

#include "assert_buffer_equals.hpp"
#include "create_execution_state.hpp"
#include "create_message.hpp"
#include "create_statement.hpp"
#include "printing.hpp"
#include "run_coroutine.hpp"
#include "test_common.hpp"
#include "test_connection.hpp"
#include "unit_netfun_maker.hpp"

using boost::mysql::blob;
using boost::mysql::client_errc;
using boost::mysql::error_code;
using boost::mysql::execution_state;
using boost::mysql::field_view;
using boost::mysql::statement;
using boost::mysql::string_view;
using boost::mysql::detail::execution_state_access;
using boost::mysql::detail::protocol_field_type;
using boost::mysql::detail::resultset_encoding;
using namespace boost::mysql::test;

namespace {

// Machinery to treat iterator and tuple overloads the same
using netfun_maker_it = netfun_maker_mem<
    void,
    test_connection,
    const statement&,
    const field_view*,
    const field_view*,
    execution_state&>;

using netfun_maker_tuple = netfun_maker_mem<
    void,
    test_connection,
    const statement&,
    const std::tuple<field_view, field_view>&,
    execution_state&>;

using start_statement_execution_fn = std::function<
    network_result<void>(test_connection&, const statement&, field_view, field_view, execution_state&)>;

start_statement_execution_fn to_common_sig(const netfun_maker_it::signature& sig)
{
    return [sig](
               test_connection& conn,
               const statement& stmt,
               field_view p1,
               field_view p2,
               execution_state& st
           ) {
        field_view params[] = {p1, p2};
        return sig(conn, stmt, std::begin(params), std::end(params), st);
    };
}

start_statement_execution_fn to_common_sig(const netfun_maker_tuple::signature& sig)
{
    return [sig](
               test_connection& conn,
               const statement& stmt,
               field_view p1,
               field_view p2,
               execution_state& st
           ) { return sig(conn, stmt, std::make_tuple(p1, p2), st); };
}

struct
{
    start_statement_execution_fn start_statement_execution;
    const char* name;
} all_fns[] = {
    {to_common_sig(netfun_maker_it::sync_errc(&test_connection::start_statement_execution)),                "sync_errc_it"},
    {to_common_sig(netfun_maker_tuple::sync_errc(&test_connection::start_statement_execution)),
     "sync_errc_tuple"                                                                                                    },
    {to_common_sig(netfun_maker_it::sync_exc(&test_connection::start_statement_execution)),                 "sync_exc_it" },
    {to_common_sig(netfun_maker_tuple::sync_exc(&test_connection::start_statement_execution)),
     "sync_exc_tuple"                                                                                                     },
    {to_common_sig(netfun_maker_it::async_errinfo(&test_connection::async_start_statement_execution)),
     "async_errinfo_it"                                                                                                   },
    {to_common_sig(netfun_maker_tuple::async_errinfo(&test_connection::async_start_statement_execution)),
     "async_errinfo_tuple"                                                                                                },
    {to_common_sig(netfun_maker_it::async_noerrinfo(&test_connection::async_start_statement_execution)),
     "async_noerrinfo_it"                                                                                                 },
    {to_common_sig(netfun_maker_tuple::async_noerrinfo(&test_connection::async_start_statement_execution)),
     "async_noerrinfo_tuple"                                                                                              },
};

// Verify that we reset the state
execution_state create_initial_state()
{
    return create_execution_state(resultset_encoding::text, {protocol_field_type::geometry}, 4);
}

BOOST_AUTO_TEST_SUITE(test_start_statement_execution)

BOOST_AUTO_TEST_CASE(success)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            execution_state st{create_initial_state()};
            auto ok_pack = create_ok_packet_message(1, 2);
            test_connection conn;
            auto stmt = create_statement(2);
            conn.stream().add_message(ok_pack);

            // Call the function
            fns.start_statement_execution(conn, stmt, field_view("test"), field_view(nullptr), st)
                .validate_no_error();

            // Verify the message we sent
            std::uint8_t expected_message[] = {
                0x015, 0x00, 0x00, 0x00, 0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00,
                0x00,  0x02, 0x01, 0xfe, 0x00, 0x06, 0x00, 0x04, 0x74, 0x65, 0x73, 0x74,
            };
            BOOST_MYSQL_ASSERT_BLOB_EQUALS(conn.stream().bytes_written(), expected_message);

            // Verify the results
            BOOST_TEST(execution_state_access::get_encoding(st) == resultset_encoding::binary);
            BOOST_TEST(st.complete());
            BOOST_TEST(execution_state_access::get_sequence_number(st) == 2u);
            BOOST_TEST(st.meta().size() == 0u);
            BOOST_TEST(st.affected_rows() == 2u);
        }
    }
}

// This covers any case where start_execution_generic would error
BOOST_AUTO_TEST_CASE(error_start_execution_generic)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            execution_state st{create_initial_state()};
            test_connection conn;
            auto stmt = create_statement(2);
            conn.stream().set_fail_count(fail_count(0, client_errc::server_unsupported));

            // Call the function
            fns.start_statement_execution(conn, stmt, field_view("test"), field_view(nullptr), st)
                .validate_error_exact(client_errc::server_unsupported);
        }
    }
}

BOOST_AUTO_TEST_CASE(error_wrong_num_params)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            execution_state st{create_initial_state()};
            test_connection conn;
            auto stmt = create_statement(3);

            // Call the function
            fns.start_statement_execution(conn, stmt, field_view("test"), field_view(nullptr), st)
                .validate_error_exact(client_errc::wrong_num_params);
        }
    }
}

// Verify that all param types we advertise work as expected
// (database_types tests use field_views)
BOOST_AUTO_TEST_CASE(tuple_parameter_types)
{
    execution_state st{create_initial_state()};
    test_connection conn;
    auto params = std::make_tuple(
        std::uint8_t(42),
        std::int16_t(-1),
        string_view("test"),
        blob({0x00, 0x01, 0x02}),
        4.2f,
        nullptr,
        boost::mysql::date(2020, 1, 2)
    );
    statement stmt{create_statement(std::tuple_size<decltype(params)>::value)};
    conn.stream().add_message(create_ok_packet_message(1));

    // Execute
    conn.start_statement_execution(stmt, params, st);

    // Verify
    std::uint8_t expected_msg[] = {
        0x3c, 0x00, 0x00, 0x00, 0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x20, 0x01,
        0x08, 0x80, 0x08, 0x00, 0xfe, 0x00, 0xfc, 0x00, 0x04, 0x00, 0x06, 0x00, 0x0a, 0x00, 0x2a, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x04, 0x74,
        0x65, 0x73, 0x74, 0x03, 0x00, 0x01, 0x02, 0x66, 0x66, 0x86, 0x40, 0x04, 0xe4, 0x07, 0x01, 0x02,
    };
    BOOST_MYSQL_ASSERT_BLOB_EQUALS(conn.stream().bytes_written(), expected_msg);
}

// Verify that we correctly perform a decay-copy of the parameters,
// relevant for deferred tokens
#ifdef BOOST_ASIO_HAS_CO_AWAIT
BOOST_AUTO_TEST_SUITE(tuple_params_copying)
struct fixture
{
    execution_state st{create_initial_state()};
    test_connection conn;
    statement stmt{create_statement(2)};

    static constexpr std::uint8_t expected_msg[]{
        0x1d, 0x00, 0x00, 0x00, 0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0xfe,
        0x00, 0x08, 0x00, 0x04, 0x74, 0x65, 0x73, 0x74, 0x2a, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    };

    fixture() { conn.stream().add_message(create_ok_packet_message(1)); }
};

BOOST_AUTO_TEST_CASE(rvalue)
{
    run_coroutine([]() -> boost::asio::awaitable<void> {
        fixture fix;

        // Deferred op
        auto aw = fix.conn.async_start_statement_execution(
            fix.stmt,
            std::make_tuple(std::string("test"), 42),
            fix.st,
            boost::asio::use_awaitable
        );
        co_await std::move(aw);

        // verify that the op had the intended effects
        BOOST_MYSQL_ASSERT_BLOB_EQUALS(fix.conn.stream().bytes_written(), fix.expected_msg);
        BOOST_TEST(fix.st.complete());
    });
}

BOOST_AUTO_TEST_CASE(lvalue)
{
    run_coroutine([]() -> boost::asio::awaitable<void> {
        fixture fix;

        // Deferred op
        auto tup = std::make_tuple(std::string("test"), 42);
        auto aw = fix.conn.async_start_statement_execution(fix.stmt, tup, fix.st, boost::asio::use_awaitable);
        tup = std::make_tuple(std::string("other"), 90);
        co_await std::move(aw);

        // verify that the op had the intended effects
        BOOST_MYSQL_ASSERT_BLOB_EQUALS(fix.conn.stream().bytes_written(), fix.expected_msg);
        BOOST_TEST(fix.st.complete());
    });
}

BOOST_AUTO_TEST_CASE(const_lvalue)
{
    run_coroutine([]() -> boost::asio::awaitable<void> {
        fixture fix;

        // Deferred op
        const auto tup = std::make_tuple(std::string("test"), 42);
        auto aw = fix.conn.async_start_statement_execution(fix.stmt, tup, fix.st, boost::asio::use_awaitable);
        co_await std::move(aw);

        // verify that the op had the intended effects
        BOOST_MYSQL_ASSERT_BLOB_EQUALS(fix.conn.stream().bytes_written(), fix.expected_msg);
        BOOST_TEST(fix.st.complete());
    });
}

BOOST_AUTO_TEST_SUITE_END()
#endif

BOOST_AUTO_TEST_SUITE_END()

}  // namespace