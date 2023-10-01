//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>

#include <boost/mysql/detail/any_stream.hpp>

#include <boost/mysql/impl/internal/channel/message_writer.hpp>
#include <boost/mysql/impl/internal/channel/write_message.hpp>

#include <boost/core/span.hpp>
#include <boost/test/unit_test.hpp>

#include <vector>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/buffer_concat.hpp"
#include "test_unit/create_frame.hpp"
#include "test_unit/test_stream.hpp"
#include "test_unit/unit_netfun_maker.hpp"

using namespace boost::mysql::detail;
using namespace boost::mysql::test;
using boost::span;
using boost::mysql::client_errc;
using boost::mysql::error_code;

BOOST_AUTO_TEST_SUITE(test_write_message)

using netfun_maker = netfun_maker_fn<void, any_stream&, message_writer&>;

struct
{
    netfun_maker::signature write;
    const char* name;
} all_fns[] = {
    {netfun_maker::sync_errc_noerrinfo(&write_message),   "sync" },
    {netfun_maker::async_noerrinfo(&async_write_message), "async"},
};

void copy(span<const std::uint8_t> from, span<std::uint8_t> to)
{
    BOOST_TEST_REQUIRE(from.size() == to.size());
    std::memcpy(to.data(), from.data(), from.size());
}

struct fixture
{
    message_writer writer{8};
    std::uint8_t seqnum{4};
    test_any_stream stream;

    test_stream& inner_stream() noexcept { return cast<test_stream>(stream); }
};

BOOST_AUTO_TEST_CASE(success)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            // Data
            fixture fix;
            const std::vector<std::uint8_t> msg{0x01, 0x02, 0x03};

            // Setup data to be written
            auto mutbuf = fix.writer.prepare_buffer(msg.size(), fix.seqnum);
            copy(msg, mutbuf);

            // Write
            fns.write(fix.stream, fix.writer).validate_no_error();

            // Verify
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(fix.inner_stream().bytes_written(), create_frame(4, msg));
            BOOST_TEST(fix.seqnum == 5);
        }
    }
}

BOOST_AUTO_TEST_CASE(empty_message)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            // Data
            fixture fix;

            // Setup data to be written
            fix.writer.prepare_buffer(0, fix.seqnum);

            // Write
            fns.write(fix.stream, fix.writer).validate_no_error();

            // Verify
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(fix.inner_stream().bytes_written(), create_empty_frame(4));
            BOOST_TEST(fix.seqnum == 5);
        }
    }
}

BOOST_AUTO_TEST_CASE(short_writes)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            // Data
            fixture fix;
            fix.inner_stream().set_write_break_size(2);  // 2 byte reads
            const std::vector<std::uint8_t> msg{0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07};

            // Setup data to be written
            auto mutbuf = fix.writer.prepare_buffer(msg.size(), fix.seqnum);
            copy(msg, mutbuf);

            // Write
            fns.write(fix.stream, fix.writer).validate_no_error();

            // Verify
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(fix.inner_stream().bytes_written(), create_frame(4, msg));
            BOOST_TEST(fix.seqnum == 5);
        }
    }
}

BOOST_AUTO_TEST_CASE(multi_frame)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            // Data
            fixture fix;
            std::vector<std::uint8_t> msg_frame_1{0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08};
            std::vector<std::uint8_t> msg_frame_2{0x11, 0x12};
            const auto msg = concat_copy(msg_frame_1, msg_frame_2);
            const auto expected_msg = concat_copy(create_frame(4, msg_frame_1), create_frame(5, msg_frame_2));

            // Setup data to be written
            auto mutbuf = fix.writer.prepare_buffer(msg.size(), fix.seqnum);
            copy(msg, mutbuf);

            // Write
            fns.write(fix.stream, fix.writer).validate_no_error();

            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(fix.inner_stream().bytes_written(), expected_msg);
            BOOST_TEST(fix.seqnum == 6);
        }
    }
}

BOOST_AUTO_TEST_CASE(error)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            // Data
            fixture fix;
            fix.inner_stream().set_fail_count(fail_count(0, client_errc::server_unsupported));
            const std::vector<std::uint8_t> msg{0x01, 0x02, 0x03};

            // Setup data to be written
            auto mutbuf = fix.writer.prepare_buffer(msg.size(), fix.seqnum);
            copy(msg, mutbuf);

            // Write
            fns.write(fix.stream, fix.writer).validate_error_exact(client_errc::server_unsupported);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
