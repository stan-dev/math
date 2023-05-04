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

#include <boost/mysql/detail/protocol/capabilities.hpp>
#include <boost/mysql/detail/protocol/db_flavor.hpp>

#include <boost/asio/buffer.hpp>
#include <boost/test/unit_test.hpp>

#include "assert_buffer_equals.hpp"
#include "create_message.hpp"
#include "test_connection.hpp"
#include "unit_netfun_maker.hpp"

using namespace boost::mysql::test;
using boost::mysql::client_errc;
using boost::mysql::common_server_errc;
using boost::mysql::error_code;
using boost::mysql::detail::capabilities;

namespace {

using netfun_maker = netfun_maker_mem<void, test_connection>;

struct
{
    netfun_maker::signature ping;
    const char* name;
} all_fns[] = {
    {netfun_maker::sync_errc(&test_connection::ping),             "sync_errc"      },
    {netfun_maker::sync_exc(&test_connection::ping),              "sync_exc"       },
    {netfun_maker::async_errinfo(&test_connection::async_ping),   "async_errinfo"  },
    {netfun_maker::async_noerrinfo(&test_connection::async_ping), "async_noerrinfo"},
};

BOOST_AUTO_TEST_SUITE(test_ping)

BOOST_AUTO_TEST_CASE(process_ping_response_)
{
    struct
    {
        const char* name;
        std::vector<std::uint8_t> message;
        error_code expected_err;
        const char* expected_msg;
    } test_cases[] = {
        {"success", create_ok_packet_body(), error_code(), ""},
        {"empty_message", {}, client_errc::incomplete_message, ""},
        {"invalid_message_type", {0xab}, client_errc::protocol_value_error, ""},
        {"bad_ok_packet", {0x00, 0x01}, client_errc::incomplete_message, ""},
        {"err_packet",
         create_err_packet_body(common_server_errc::er_bad_db_error, "abc"),
         common_server_errc::er_bad_db_error,
         "abc"},
        {"bad_err_packet", {0xff, 0x01}, client_errc::incomplete_message, ""},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            boost::mysql::diagnostics diag;
            auto err = boost::mysql::detail::process_ping_response(
                boost::asio::buffer(tc.message),
                capabilities(),
                boost::mysql::detail::db_flavor::mariadb,
                diag
            );

            BOOST_TEST(err == tc.expected_err);
            BOOST_TEST(diag.server_message() == tc.expected_msg);
        }
    }
}

BOOST_AUTO_TEST_CASE(success)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            test_connection conn;
            conn.stream().add_message(create_ok_packet_message(1, 2, 3, 4, 5));

            // Call the function
            fns.ping(conn).validate_no_error();

            // Verify the message we sent
            std::uint8_t expected_message[] = {0x01, 0x00, 0x00, 0x00, 0x0e};
            BOOST_MYSQL_ASSERT_BLOB_EQUALS(conn.stream().bytes_written(), expected_message);
        }
    }
}

BOOST_AUTO_TEST_CASE(error_network)
{
    for (auto fns : all_fns)
    {
        for (int i = 0; i <= 1; ++i)
        {
            BOOST_TEST_CONTEXT(fns.name << " in network transfer " << i)
            {
                test_connection conn;
                conn.stream().set_fail_count(fail_count(i, common_server_errc::er_aborting_connection));

                // Call the function
                fns.ping(conn).validate_error_exact(common_server_errc::er_aborting_connection);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(error_response)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            test_connection conn;
            conn.stream().add_message(
                create_err_packet_message(1, common_server_errc::er_bad_db_error, "my_message")
            );

            // Call the function
            fns.ping(conn).validate_error_exact(common_server_errc::er_bad_db_error, "my_message");
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace