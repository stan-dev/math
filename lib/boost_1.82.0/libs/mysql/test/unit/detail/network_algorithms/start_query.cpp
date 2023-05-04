//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/common_server_errc.hpp>

#include <boost/mysql/detail/auxiliar/access_fwd.hpp>

#include <boost/test/unit_test.hpp>

#include "assert_buffer_equals.hpp"
#include "create_execution_state.hpp"
#include "create_message.hpp"
#include "printing.hpp"
#include "test_common.hpp"
#include "test_connection.hpp"
#include "unit_netfun_maker.hpp"

using boost::mysql::blob;
using boost::mysql::common_server_errc;
using boost::mysql::error_code;
using boost::mysql::execution_state;
using boost::mysql::string_view;
using boost::mysql::detail::execution_state_access;
using boost::mysql::detail::protocol_field_type;
using boost::mysql::detail::resultset_encoding;
using namespace boost::mysql::test;

namespace {

using netfun_maker = netfun_maker_mem<void, test_connection, string_view, execution_state&>;

struct
{
    netfun_maker::signature start_query;
    const char* name;
} all_fns[] = {
    {netfun_maker::sync_errc(&test_connection::start_query),             "sync_errc"      },
    {netfun_maker::sync_exc(&test_connection::start_query),              "sync_exc"       },
    {netfun_maker::async_errinfo(&test_connection::async_start_query),   "async_errinfo"  },
    {netfun_maker::async_noerrinfo(&test_connection::async_start_query), "async_noerrinfo"},
};

// Verify that we reset the state
execution_state create_initial_state()
{
    return create_execution_state(resultset_encoding::binary, {protocol_field_type::geometry}, 4);
}

BOOST_AUTO_TEST_SUITE(test_start_query)

BOOST_AUTO_TEST_CASE(success)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            execution_state st{create_initial_state()};
            test_connection conn;
            conn.stream().add_message(create_ok_packet_message(1, 2));

            // Call the function
            fns.start_query(conn, "SELECT 1", st).validate_no_error();

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
            BOOST_TEST(execution_state_access::get_encoding(st) == resultset_encoding::text);
            BOOST_TEST(st.complete());
            BOOST_TEST(execution_state_access::get_sequence_number(st) == 2u);
            BOOST_TEST(st.meta().size() == 0u);
            BOOST_TEST(st.affected_rows() == 2u);
        }
    }
}

BOOST_AUTO_TEST_CASE(error)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            execution_state st{create_initial_state()};
            test_connection conn;
            conn.stream().set_fail_count(fail_count(0, common_server_errc::er_aborting_connection));

            // Call the function
            fns.start_query(conn, "SELECT 1", st)
                .validate_error_exact(common_server_errc::er_aborting_connection);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace