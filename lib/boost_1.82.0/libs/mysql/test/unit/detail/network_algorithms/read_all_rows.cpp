//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/row.hpp>

#include <boost/mysql/detail/auxiliar/access_fwd.hpp>
#include <boost/mysql/detail/network_algorithms/read_all_rows.hpp>
#include <boost/mysql/detail/protocol/resultset_encoding.hpp>

#include <boost/asio/buffer.hpp>
#include <boost/test/unit_test.hpp>

#include "create_execution_state.hpp"
#include "create_message.hpp"
#include "test_channel.hpp"
#include "test_common.hpp"
#include "unit_netfun_maker.hpp"

using boost::mysql::client_errc;
using boost::mysql::error_code;
using boost::mysql::execution_state;
using boost::mysql::rows;
using boost::mysql::detail::async_read_all_rows;
using boost::mysql::detail::execution_state_access;
using boost::mysql::detail::protocol_field_type;
using boost::mysql::detail::read_all_rows;
using boost::mysql::detail::resultset_encoding;
using namespace boost::mysql::test;

namespace {

using netfun_maker = netfun_maker_fn<void, test_channel&, execution_state&, rows&>;

struct
{
    typename netfun_maker::signature read_all_rows;
    const char* name;
} all_fns[] = {
    {netfun_maker::sync_errc(&read_all_rows<test_stream>),           "sync" },
    {netfun_maker::async_errinfo(&async_read_all_rows<test_stream>), "async"}
};

// Verify that we clear any previous result
rows make_initial_rows() { return makerows(2, 42, nullptr, 4.5f, "abc"); }

BOOST_AUTO_TEST_SUITE(test_read_all_rows)

BOOST_AUTO_TEST_CASE(success_row_row_eof)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            auto row1 = create_message(4, {0x00, 0x00, 0x03, 0x6d, 0x69, 0x6e, 0x6d, 0x07});
            auto row2 = create_message(5, {0x00, 0x08, 0x03, 0x6d, 0x61, 0x78});
            auto ok_packet = create_message(6, {0xfe, 0x01, 0x06, 0x02, 0x00, 0x09, 0x00, 0x02, 0x61, 0x62});
            auto st = create_execution_state(
                resultset_encoding::binary,
                {protocol_field_type::var_string, protocol_field_type::short_},
                4  // seqnum
            );
            auto chan = create_channel(concat_copy(row1, row2, ok_packet), 1024);
            chan.shared_fields().emplace_back("abc");  // from previous call
            rows rws = make_initial_rows();

            fns.read_all_rows(chan, st, rws).validate_no_error();
            BOOST_TEST_REQUIRE(rws.size() == 2u);
            BOOST_TEST(rws[0] == makerow("min", 1901));
            BOOST_TEST(rws[1] == makerow("max", nullptr));
            BOOST_TEST(st.complete());
            BOOST_TEST(st.affected_rows() == 1u);
            BOOST_TEST(st.last_insert_id() == 6u);
            BOOST_TEST(st.warning_count() == 9u);
            BOOST_TEST(st.info() == "ab");
            BOOST_TEST(chan.shared_sequence_number() == 0u);  // not used
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
            auto ok_packet = create_message(6, {0xfe, 0x01, 0x06, 0x02, 0x00, 0x09, 0x00, 0x02, 0x61, 0x62});
            auto st = create_execution_state(
                resultset_encoding::binary,
                {protocol_field_type::var_string, protocol_field_type::short_},
                4  // seqnum
            );
            auto chan = create_channel({}, 1024);
            chan.lowest_layer().add_message(row1);
            chan.lowest_layer().add_message(concat_copy(row2, ok_packet));
            chan.shared_fields().emplace_back("abc");  // from previous call
            rows rws = make_initial_rows();

            fns.read_all_rows(chan, st, rws).validate_no_error();
            BOOST_TEST_REQUIRE(rws.size() == 2u);
            BOOST_TEST(rws[0] == makerow("min", 1901));
            BOOST_TEST(rws[1] == makerow("max", nullptr));
            BOOST_TEST(st.complete());
            BOOST_TEST(st.affected_rows() == 1u);
            BOOST_TEST(st.last_insert_id() == 6u);
            BOOST_TEST(st.warning_count() == 9u);
            BOOST_TEST(st.info() == "ab");
            BOOST_TEST(chan.shared_sequence_number() == 0u);  // not used
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
            auto ok_packet = create_message(5, {0xfe, 0x01, 0x06, 0x02, 0x00, 0x09, 0x00, 0x02, 0x61, 0x62});
            auto st = create_execution_state(
                resultset_encoding::binary,
                {protocol_field_type::var_string, protocol_field_type::short_},
                4  // seqnum
            );
            auto chan = create_channel({}, 1024);
            chan.lowest_layer().add_message(row1);
            chan.lowest_layer().add_message(ok_packet);
            chan.shared_fields().emplace_back("abc");  // from previous call
            rows rws = make_initial_rows();

            // row
            fns.read_all_rows(chan, st, rws).validate_no_error();
            BOOST_TEST_REQUIRE(rws.size() == 1u);
            BOOST_TEST(rws[0] == makerow("min", 1901));
            BOOST_TEST(st.complete());
            BOOST_TEST(st.affected_rows() == 1u);
            BOOST_TEST(st.last_insert_id() == 6u);
            BOOST_TEST(st.warning_count() == 9u);
            BOOST_TEST(st.info() == "ab");
            BOOST_TEST(chan.shared_sequence_number() == 0u);  // not used
        }
    }
}

BOOST_AUTO_TEST_CASE(success_eof)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            auto ok_packet = create_message(4, {0xfe, 0x01, 0x06, 0x02, 0x00, 0x09, 0x00, 0x02, 0x61, 0x62});
            auto st = create_execution_state(
                resultset_encoding::binary,
                {protocol_field_type::var_string, protocol_field_type::short_},
                4  // seqnum
            );
            auto chan = create_channel(ok_packet, 1024);
            chan.shared_fields().emplace_back("abc");  // from previous call
            rows rws = make_initial_rows();

            fns.read_all_rows(chan, st, rws).validate_no_error();
            BOOST_TEST(rws.size() == 0u);
            BOOST_TEST(st.complete());
            BOOST_TEST(st.affected_rows() == 1u);
            BOOST_TEST(st.last_insert_id() == 6u);
            BOOST_TEST(st.warning_count() == 9u);
            BOOST_TEST(st.info() == "ab");
            BOOST_TEST(chan.shared_sequence_number() == 0u);  // not used
        }
    }
}

// caught as failing by an integ test
BOOST_AUTO_TEST_CASE(success_eof_shared_fields_empty)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            auto ok_packet = create_message(4, {0xfe, 0x01, 0x06, 0x02, 0x00, 0x09, 0x00, 0x02, 0x61, 0x62});
            auto st = create_execution_state(
                resultset_encoding::binary,
                {protocol_field_type::var_string, protocol_field_type::short_},
                4  // seqnum
            );
            auto chan = create_channel(ok_packet, 1024);
            rows rws = make_initial_rows();

            fns.read_all_rows(chan, st, rws).validate_no_error();
            BOOST_TEST(rws.size() == 0u);
            BOOST_TEST(st.complete());
            BOOST_TEST(st.affected_rows() == 1u);
            BOOST_TEST(st.last_insert_id() == 6u);
            BOOST_TEST(st.warning_count() == 9u);
            BOOST_TEST(st.info() == "ab");
            BOOST_TEST(chan.shared_sequence_number() == 0u);  // not used
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
            test_channel chan = create_channel();
            rows rws = make_initial_rows();

            fns.read_all_rows(chan, st, rws).validate_no_error();
            BOOST_TEST(rws.empty());
            BOOST_TEST(st.complete());

            // Doing it again works, too
            fns.read_all_rows(chan, st, rws).validate_no_error();
            BOOST_TEST(rws.empty());
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
            test_channel chan = create_channel();
            rows rws = make_initial_rows();
            chan.lowest_layer().set_fail_count(fail_count(0, client_errc::server_unsupported));

            fns.read_all_rows(chan, st, rws).validate_error_exact(client_errc::server_unsupported);
            BOOST_TEST(rws.empty());
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
            test_channel chan = create_channel();
            rows rws = make_initial_rows();
            chan.lowest_layer().add_message(r);

            // deserialize row error
            fns.read_all_rows(chan, st, rws).validate_error_exact(client_errc::incomplete_message);
            BOOST_TEST(rws.empty());
            BOOST_TEST(!st.complete());
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace