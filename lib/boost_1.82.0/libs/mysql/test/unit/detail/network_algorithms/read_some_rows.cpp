//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/row_view.hpp>

#include <boost/mysql/detail/auxiliar/access_fwd.hpp>
#include <boost/mysql/detail/protocol/resultset_encoding.hpp>

#include <boost/asio/buffer.hpp>
#include <boost/test/unit_test.hpp>

#include "create_execution_state.hpp"
#include "create_message.hpp"
#include "test_channel.hpp"
#include "test_common.hpp"
#include "test_connection.hpp"
#include "unit_netfun_maker.hpp"

using boost::mysql::client_errc;
using boost::mysql::error_code;
using boost::mysql::execution_state;
using boost::mysql::rows_view;
using boost::mysql::detail::execution_state_access;
using boost::mysql::detail::protocol_field_type;
using boost::mysql::detail::resultset_encoding;
using namespace boost::mysql::test;

namespace {

test_channel& get_channel(test_connection& conn) noexcept
{
    return boost::mysql::detail::connection_access::get_channel(conn);
}

using netfun_maker = netfun_maker_mem<rows_view, test_connection, execution_state&>;

struct
{
    typename netfun_maker::signature read_some_rows;
    const char* name;
} all_fns[] = {
    {netfun_maker::sync_errc(&test_connection::read_some_rows),           "sync" },
    {netfun_maker::async_errinfo(&test_connection::async_read_some_rows), "async"},
};

BOOST_AUTO_TEST_SUITE(test_read_some_rows)

BOOST_AUTO_TEST_CASE(success_row_row_eof)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            auto row1 = create_message(4, {0x00, 0x00, 0x03, 0x6d, 0x69, 0x6e, 0x6d, 0x07});
            auto row2 = create_message(5, {0x00, 0x08, 0x03, 0x6d, 0x61, 0x78});
            auto ok_packet = create_eof_packet_message(6, 1, 6, 0, 9, "ab");
            auto st = create_execution_state(
                resultset_encoding::binary,
                {protocol_field_type::var_string, protocol_field_type::short_},
                4  // seqnum
            );
            test_connection conn;
            get_channel(conn).shared_fields().emplace_back("abc");  // from previous call
            conn.stream().add_message(concat_copy(row1, row2, ok_packet));

            rows_view rv = fns.read_some_rows(conn, st).get();
            BOOST_TEST_REQUIRE(rv.size() == 2u);
            BOOST_TEST(rv[0] == makerow("min", 1901));
            BOOST_TEST(rv[1] == makerow("max", nullptr));
            BOOST_TEST(st.complete());
            BOOST_TEST(st.affected_rows() == 1u);
            BOOST_TEST(st.last_insert_id() == 6u);
            BOOST_TEST(st.warning_count() == 9u);
            BOOST_TEST(st.info() == "ab");
            BOOST_TEST(get_channel(conn).shared_sequence_number() == 0u);  // not used
        }
    }
}

BOOST_AUTO_TEST_CASE(success_row_row_eof_separate)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            auto row1 = create_message(4, {0x00, 0x00, 0x03, 0x6d, 0x69, 0x6e, 0x6d, 0x07});
            auto row2 = create_message(5, {0x00, 0x08, 0x03, 0x6d, 0x61, 0x78});
            auto ok_packet = create_eof_packet_message(6, 1, 6, 0, 9, "ab");
            auto st = create_execution_state(
                resultset_encoding::binary,
                {protocol_field_type::var_string, protocol_field_type::short_},
                4  // seqnum
            );
            test_connection conn;
            get_channel(conn).shared_fields().emplace_back("abc");  // from previous call
            conn.stream().add_message(row1);
            conn.stream().add_message(concat_copy(row2, ok_packet));

            // 1st read
            rows_view rv = fns.read_some_rows(conn, st).get();
            BOOST_TEST_REQUIRE(rv.size() == 1u);
            BOOST_TEST(rv[0] == makerow("min", 1901));
            BOOST_TEST(!st.complete());

            // 2nd read
            rv = fns.read_some_rows(conn, st).get();
            BOOST_TEST_REQUIRE(rv.size() == 1u);
            BOOST_TEST(rv[0] == makerow("max", nullptr));
            BOOST_TEST(st.complete());
            BOOST_TEST(st.affected_rows() == 1u);
            BOOST_TEST(st.last_insert_id() == 6u);
            BOOST_TEST(st.warning_count() == 9u);
            BOOST_TEST(st.info() == "ab");
            BOOST_TEST(get_channel(conn).shared_sequence_number() == 0u);  // not used
        }
    }
}

BOOST_AUTO_TEST_CASE(success_row_eof_separate)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            auto row1 = create_message(4, {0x00, 0x00, 0x03, 0x6d, 0x69, 0x6e, 0x6d, 0x07});
            auto ok_packet = create_eof_packet_message(5, 1, 6, 0, 9, "ab");
            auto st = create_execution_state(
                resultset_encoding::binary,
                {protocol_field_type::var_string, protocol_field_type::short_},
                4  // seqnum
            );
            test_connection conn;
            get_channel(conn).shared_fields().emplace_back("abc");  // from previous call
            conn.stream().add_message(row1);
            conn.stream().add_message(ok_packet);

            // row
            rows_view rv = fns.read_some_rows(conn, st).get();
            BOOST_TEST_REQUIRE(rv.size() == 1u);
            BOOST_TEST(rv[0] == makerow("min", 1901));
            BOOST_TEST(!st.complete());

            // eof
            rv = fns.read_some_rows(conn, st).get();
            BOOST_TEST(rv.size() == 0u);
            BOOST_TEST(st.complete());
            BOOST_TEST(st.affected_rows() == 1u);
            BOOST_TEST(st.last_insert_id() == 6u);
            BOOST_TEST(st.warning_count() == 9u);
            BOOST_TEST(st.info() == "ab");
            BOOST_TEST(get_channel(conn).shared_sequence_number() == 0u);  // not used
        }
    }
}

BOOST_AUTO_TEST_CASE(success_eof)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            auto ok_packet = create_eof_packet_message(4, 1, 6, 0, 9, "ab");
            auto st = create_execution_state(
                resultset_encoding::binary,
                {protocol_field_type::var_string, protocol_field_type::short_},
                4  // seqnum
            );
            test_connection conn;
            get_channel(conn).shared_fields().emplace_back("abc");  // from previous call
            conn.stream().add_message(ok_packet);

            rows_view rv = fns.read_some_rows(conn, st).get();
            BOOST_TEST(rv.size() == 0u);
            BOOST_TEST(st.complete());
            BOOST_TEST(st.affected_rows() == 1u);
            BOOST_TEST(st.last_insert_id() == 6u);
            BOOST_TEST(st.warning_count() == 9u);
            BOOST_TEST(st.info() == "ab");
            BOOST_TEST(get_channel(conn).shared_sequence_number() == 0u);  // not used
        }
    }
}

BOOST_AUTO_TEST_CASE(resultset_already_complete)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            auto st = create_execution_state(resultset_encoding::text, {});
            execution_state_access::complete(st, boost::mysql::detail::ok_packet{});
            test_connection conn;

            rows_view rv = fns.read_some_rows(conn, st).get();
            BOOST_TEST(rv.empty());
            BOOST_TEST(st.complete());

            // Doing it again works, too
            rv = fns.read_some_rows(conn, st).get();
            BOOST_TEST(rv.empty());
            BOOST_TEST(st.complete());
        }
    }
}

BOOST_AUTO_TEST_CASE(error_reading_row)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            auto st = create_execution_state(resultset_encoding::text, {});
            test_connection conn;
            conn.stream().set_fail_count(fail_count(0, client_errc::server_unsupported));

            fns.read_some_rows(conn, st).validate_error_exact(client_errc::server_unsupported);
            BOOST_TEST(!st.complete());
        }
    }
}

BOOST_AUTO_TEST_CASE(error_deserializing_row)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            auto r = create_message(0, {0x00});  // invalid row
            auto st = create_execution_state(resultset_encoding::binary, {protocol_field_type::var_string});
            test_connection conn;
            conn.stream().add_message(r);

            // deserialize row error
            fns.read_some_rows(conn, st).validate_error_exact(client_errc::incomplete_message);
            BOOST_TEST(!st.complete());
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace