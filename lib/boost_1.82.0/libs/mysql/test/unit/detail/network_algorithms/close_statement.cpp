//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/statement.hpp>

#include <boost/asio/awaitable.hpp>
#include <boost/test/unit_test.hpp>

#include "assert_buffer_equals.hpp"
#include "create_statement.hpp"
#include "run_coroutine.hpp"
#include "test_connection.hpp"
#include "unit_netfun_maker.hpp"

using namespace boost::mysql::test;
using boost::mysql::common_server_errc;
using boost::mysql::error_code;
using boost::mysql::statement;

namespace {

using netfun_maker = netfun_maker_mem<void, test_connection, const statement&>;

struct
{
    netfun_maker::signature close_statement;
    const char* name;
} all_fns[] = {
    {netfun_maker::sync_errc(&test_connection::close_statement),             "sync_errc"      },
    {netfun_maker::sync_exc(&test_connection::close_statement),              "sync_exc"       },
    {netfun_maker::async_errinfo(&test_connection::async_close_statement),   "async_errinfo"  },
    {netfun_maker::async_noerrinfo(&test_connection::async_close_statement), "async_noerrinfo"},
};

BOOST_AUTO_TEST_SUITE(test_close_statement)

constexpr std::uint8_t expected_message[]{0x05, 0x00, 0x00, 0x00, 0x19, 0x03, 0x00, 0x00, 0x00};

struct fixture
{
    test_connection conn;
    statement stmt{create_statement(2, 3)};
};

BOOST_AUTO_TEST_CASE(success)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;

            // Call the function
            fns.close_statement(fix.conn, fix.stmt).validate_no_error();

            // Verify the message we sent
            BOOST_MYSQL_ASSERT_BLOB_EQUALS(fix.conn.stream().bytes_written(), expected_message);
        }
    }
}

BOOST_AUTO_TEST_CASE(error_network)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.conn.stream().set_fail_count(fail_count(0, common_server_errc::er_aborting_connection));

            // Call the function
            fns.close_statement(fix.conn, fix.stmt)
                .validate_error_exact(common_server_errc::er_aborting_connection);
        }
    }
}

// Verify that we don't require the passed-in statement to be alive. Only
// relevant for deferred tokens.
#ifdef BOOST_ASIO_HAS_CO_AWAIT
BOOST_AUTO_TEST_CASE(statement_handle)
{
    run_coroutine([]() -> boost::asio::awaitable<void> {
        fixture fix;

        // Deferred op
        auto aw = fix.conn.async_close_statement(statement(fix.stmt), boost::asio::use_awaitable);
        co_await std::move(aw);

        // verify that the op had the intended effects
        BOOST_MYSQL_ASSERT_BLOB_EQUALS(fix.conn.stream().bytes_written(), expected_message);
    });
}
#endif

BOOST_AUTO_TEST_SUITE_END()

}  // namespace