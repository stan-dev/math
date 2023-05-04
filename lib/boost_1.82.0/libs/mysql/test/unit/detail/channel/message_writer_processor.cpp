//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/error_code.hpp>

#include <boost/mysql/detail/channel/message_writer_processor.hpp>
#include <boost/mysql/detail/protocol/capabilities.hpp>
#include <boost/mysql/detail/protocol/common_messages.hpp>
#include <boost/mysql/detail/protocol/deserialization_context.hpp>
#include <boost/mysql/detail/protocol/serialization.hpp>

#include <boost/asio/buffer.hpp>
#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <cstdint>
#include <vector>

#include "assert_buffer_equals.hpp"
#include "buffer_concat.hpp"
#include "create_message.hpp"

using boost::asio::buffer;
using boost::mysql::error_code;
using boost::mysql::detail::message_writer_processor;
using boost::mysql::test::concat_copy;

namespace {

void check_header(
    boost::asio::const_buffer buff,
    std::uint8_t expected_seqnum,
    std::size_t expected_size
)
{
    BOOST_TEST_REQUIRE(buff.size() == 4u);
    boost::mysql::detail::deserialization_context ctx(buff, boost::mysql::detail::capabilities());
    boost::mysql::detail::packet_header header;
    auto err = boost::mysql::detail::deserialize_message(ctx, header);
    BOOST_TEST(err == error_code());
    BOOST_TEST(header.sequence_number == expected_seqnum);
    BOOST_TEST(header.packet_size.value == expected_size);
}

BOOST_AUTO_TEST_SUITE(test_message_writer_processor)

BOOST_AUTO_TEST_CASE(regular_message)
{
    message_writer_processor processor(8);
    std::vector<std::uint8_t> msg_body{0x01, 0x02, 0x03};
    std::uint8_t seqnum = 2;

    // Operation start
    processor.reset(buffer(msg_body), seqnum);
    BOOST_TEST(!processor.is_complete());

    // Prepare chunk
    auto buffers = processor.prepare_next_chunk();
    check_header(buffers[0], 2, 3);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buffers[1], buffer(msg_body));

    // On write successful
    processor.on_bytes_written();
    BOOST_TEST(seqnum == 3u);
    BOOST_TEST(processor.is_complete());
}

BOOST_AUTO_TEST_CASE(empty_message)
{
    message_writer_processor processor(8);
    std::vector<std::uint8_t> msg_body{};
    std::uint8_t seqnum = 2;

    // Operation start
    processor.reset(buffer(msg_body), seqnum);
    BOOST_TEST(!processor.is_complete());

    // Prepare chunk
    auto buffers = processor.prepare_next_chunk();
    check_header(buffers[0], 2, 0);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buffers[1], buffer(msg_body));

    // On write successful
    processor.on_bytes_written();
    BOOST_TEST(seqnum == 3u);
    BOOST_TEST(processor.is_complete());
}

BOOST_AUTO_TEST_CASE(message_with_max_frame_size_length)
{
    message_writer_processor processor(8);
    std::vector<std::uint8_t> msg_body{0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08};
    std::uint8_t seqnum = 2;

    // Operation start
    processor.reset(buffer(msg_body), seqnum);
    BOOST_TEST(!processor.is_complete());

    // Prepare chunk
    auto buffers = processor.prepare_next_chunk();
    check_header(buffers[0], 2, 8);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buffers[1], buffer(msg_body));

    // On write successful
    processor.on_bytes_written();
    BOOST_TEST(!processor.is_complete());

    // Prepare next chunk
    buffers = processor.prepare_next_chunk();
    check_header(buffers[0], 3, 0);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buffers[1], boost::asio::const_buffer());

    // On write successful
    processor.on_bytes_written();
    BOOST_TEST(seqnum == 4u);
    BOOST_TEST(processor.is_complete());
}

BOOST_AUTO_TEST_CASE(multiframe_message)
{
    message_writer_processor processor(8);
    std::vector<std::uint8_t> msg_frame_1{0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08};
    std::vector<std::uint8_t> msg_frame_2{0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18};
    std::vector<std::uint8_t> msg_frame_3{0x21};
    auto msg = concat_copy(msg_frame_1, msg_frame_2, msg_frame_3);
    std::uint8_t seqnum = 2;

    // Operation start
    processor.reset(buffer(msg), seqnum);
    BOOST_TEST(!processor.is_complete());

    // Prepare chunk 1
    auto buffers = processor.prepare_next_chunk();
    check_header(buffers[0], 2, 8);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buffers[1], buffer(msg_frame_1));

    // On write successful 1
    processor.on_bytes_written();
    BOOST_TEST(!processor.is_complete());

    // Prepare next chunk 2
    buffers = processor.prepare_next_chunk();
    check_header(buffers[0], 3, 8);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buffers[1], buffer(msg_frame_2));

    // On write successful 2
    processor.on_bytes_written();
    BOOST_TEST(!processor.is_complete());

    // Prepare next chunk 3
    buffers = processor.prepare_next_chunk();
    check_header(buffers[0], 4, 1);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buffers[1], buffer(msg_frame_3));

    // On write successful
    processor.on_bytes_written();
    BOOST_TEST(seqnum == 5u);
    BOOST_TEST(processor.is_complete());
}

BOOST_AUTO_TEST_CASE(multiframe_message_with_max_frame_size)
{
    message_writer_processor processor(8);
    std::vector<std::uint8_t> msg_frame_1{0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08};
    std::vector<std::uint8_t> msg_frame_2{0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18};
    auto msg = concat_copy(msg_frame_1, msg_frame_2);
    std::uint8_t seqnum = 2;

    // Operation start
    processor.reset(buffer(msg), seqnum);
    BOOST_TEST(!processor.is_complete());

    // Prepare chunk 1
    auto buffers = processor.prepare_next_chunk();
    check_header(buffers[0], 2, 8);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buffers[1], buffer(msg_frame_1));

    // On write successful 1
    processor.on_bytes_written();
    BOOST_TEST(!processor.is_complete());

    // Prepare next chunk 2
    buffers = processor.prepare_next_chunk();
    check_header(buffers[0], 3, 8);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buffers[1], buffer(msg_frame_2));

    // On write successful 2
    processor.on_bytes_written();
    BOOST_TEST(!processor.is_complete());

    // Prepare next chunk 3
    buffers = processor.prepare_next_chunk();
    check_header(buffers[0], 4, 0);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buffers[1], boost::asio::const_buffer());

    // On write successful
    processor.on_bytes_written();
    BOOST_TEST(seqnum == 5u);
    BOOST_TEST(processor.is_complete());
}

BOOST_AUTO_TEST_CASE(seqnum_overflow)
{
    message_writer_processor processor(8);
    std::vector<std::uint8_t> msg{0x01, 0x02};
    std::uint8_t seqnum = 0xff;

    // Operation start
    processor.reset(buffer(msg), seqnum);
    BOOST_TEST(!processor.is_complete());

    // Prepare chunk
    auto buffers = processor.prepare_next_chunk();
    check_header(buffers[0], 0xff, 2);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buffers[1], buffer(msg));

    // On write successful
    processor.on_bytes_written();
    BOOST_TEST(seqnum == 0);
    BOOST_TEST(processor.is_complete());
}

BOOST_AUTO_TEST_CASE(several_messages)
{
    message_writer_processor processor(8);
    std::vector<std::uint8_t> msg_1{0x01, 0x02};
    std::vector<std::uint8_t> msg_2{0x04, 0x05, 0x06};
    std::uint8_t seqnum_1 = 2;
    std::uint8_t seqnum_2 = 42;

    // Operation start for message 1
    processor.reset(buffer(msg_1), seqnum_1);
    BOOST_TEST(!processor.is_complete());

    // Prepare chunk for message 1
    auto buffers = processor.prepare_next_chunk();
    check_header(buffers[0], 2, 2);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buffers[1], buffer(msg_1));

    // On write successful message 1
    processor.on_bytes_written();
    BOOST_TEST(seqnum_1 == 3u);
    BOOST_TEST(processor.is_complete());

    // Operation start for message 2
    processor.reset(buffer(msg_2), seqnum_2);
    BOOST_TEST(!processor.is_complete());

    // Prepare chunk for message 2
    buffers = processor.prepare_next_chunk();
    check_header(buffers[0], 42, 3);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buffers[1], buffer(msg_2));

    // On write successful message 2
    processor.on_bytes_written();
    BOOST_TEST(seqnum_2 == 43u);
    BOOST_TEST(processor.is_complete());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace
