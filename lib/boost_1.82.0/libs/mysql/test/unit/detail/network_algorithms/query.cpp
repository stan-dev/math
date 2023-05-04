//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/column_type.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/results.hpp>

#include <boost/mysql/detail/auxiliar/access_fwd.hpp>
#include <boost/mysql/detail/protocol/constants.hpp>
#include <boost/mysql/detail/protocol/resultset_encoding.hpp>

#include <boost/test/unit_test.hpp>

#include "assert_buffer_equals.hpp"
#include "create_execution_state.hpp"
#include "create_message.hpp"
#include "printing.hpp"
#include "run_coroutine.hpp"
#include "test_common.hpp"
#include "test_connection.hpp"
#include "unit_netfun_maker.hpp"

using boost::mysql::blob;
using boost::mysql::column_type;
using boost::mysql::common_server_errc;
using boost::mysql::error_code;
using boost::mysql::results;
using boost::mysql::string_view;
using boost::mysql::detail::connection_access;
using boost::mysql::detail::protocol_field_type;
using boost::mysql::detail::results_access;
using boost::mysql::detail::resultset_encoding;
using namespace boost::mysql::test;

namespace {

using netfun_maker = netfun_maker_mem<void, test_connection, string_view, results&>;

struct
{
    netfun_maker::signature query;
    const char* name;
} all_fns[] = {
    {netfun_maker::sync_errc(&test_connection::query),             "sync_errc"      },
    {netfun_maker::sync_exc(&test_connection::query),              "sync_exc"       },
    {netfun_maker::async_errinfo(&test_connection::async_query),   "async_errinfo"  },
    {netfun_maker::async_noerrinfo(&test_connection::async_query), "async_noerrinfo"},
};

// Verify that we reset the results object
results create_initial_results()
{
    results res;
    results_access::get_rows(res) = makerows(1, 42, "abc");
    results_access::get_state(res
    ) = create_execution_state(resultset_encoding::binary, {protocol_field_type::geometry}, 4);
    return res;
}

BOOST_AUTO_TEST_SUITE(test_query)

BOOST_AUTO_TEST_CASE(success)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            auto result = create_initial_results();
            test_connection conn;
            conn.stream().add_message(create_ok_packet_message(1, 2, 3, 4, 5, "info"));

            // Call the function
            fns.query(conn, "SELECT 1", result).validate_no_error();

            // Verify the message we sent
            std::uint8_t expected_message[] = {
                0x09,
                0x00,
                0x00,
                0x00,
                0x03,
                0x53,
                0x45,
                0x4c,
                0x45,
                0x43,
                0x54,
                0x20,
                0x31,
            };
            BOOST_MYSQL_ASSERT_BLOB_EQUALS(conn.stream().bytes_written(), expected_message);

            // Verify the results
            BOOST_TEST(result.meta().size() == 0u);
            BOOST_TEST(result.affected_rows() == 2u);
            BOOST_TEST(result.last_insert_id() == 3u);
            BOOST_TEST(result.warning_count() == 5u);
            BOOST_TEST(result.info() == "info");
        }
    }
}

BOOST_AUTO_TEST_CASE(error_start_query)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            auto result = create_initial_results();
            test_connection conn;
            conn.stream().set_fail_count(fail_count(0, common_server_errc::er_aborting_connection));

            // Call the function
            fns.query(conn, "SELECT 1", result)
                .validate_error_exact(common_server_errc::er_aborting_connection);
        }
    }
}

BOOST_AUTO_TEST_CASE(error_read_all_rows)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            auto result = create_initial_results();
            test_connection conn;
            conn.stream().add_message(create_message(1, {0x01}));  // Response OK, 1 metadata packet
            conn.stream().add_message(create_coldef_message(2, protocol_field_type::geometry));
            conn.stream().set_fail_count(fail_count(4, common_server_errc::er_aborting_connection));

            // Call the function
            fns.query(conn, "SELECT 1", result)
                .validate_error_exact(common_server_errc::er_aborting_connection);

            // Ensure we successfully ran the start_query
            BOOST_TEST_REQUIRE(results_access::get_state(result).meta().size() == 1u);
            BOOST_TEST(results_access::get_state(result).meta()[0].type() == column_type::geometry);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace