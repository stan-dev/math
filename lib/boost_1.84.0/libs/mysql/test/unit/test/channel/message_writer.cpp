//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/impl/internal/channel/message_writer.hpp>

#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <cstdint>
#include <vector>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/buffer_concat.hpp"
#include "test_unit/create_frame.hpp"

using namespace boost::mysql::detail;
using namespace boost::mysql::test;
using boost::span;

namespace {

BOOST_AUTO_TEST_SUITE(test_message_writer)

void copy(span<const std::uint8_t> from, span<std::uint8_t> to)
{
    BOOST_TEST_REQUIRE(from.size() == to.size());
    std::memcpy(to.data(), from.data(), from.size());
}

BOOST_AUTO_TEST_SUITE(chunk_processor_)
BOOST_AUTO_TEST_CASE(reset)
{
    chunk_processor chunk_proc;
    BOOST_TEST(chunk_proc.done());
}

BOOST_AUTO_TEST_CASE(zero_to_size_steps)
{
    chunk_processor chunk_proc;
    std::vector<std::uint8_t> buff(10, 0);

    chunk_proc.reset(0, 10);
    BOOST_TEST(!chunk_proc.done());
    auto chunk = chunk_proc.get_chunk(buff);
    BOOST_TEST(chunk.data() == buff.data());
    BOOST_TEST(chunk.size() == 10u);

    chunk_proc.on_bytes_written(3);
    BOOST_TEST(!chunk_proc.done());
    chunk = chunk_proc.get_chunk(buff);
    BOOST_TEST(chunk.data() == buff.data() + 3);
    BOOST_TEST(chunk.size() == 7u);

    chunk_proc.on_bytes_written(6);
    BOOST_TEST(!chunk_proc.done());
    chunk = chunk_proc.get_chunk(buff);
    BOOST_TEST(chunk.data() == buff.data() + 9);
    BOOST_TEST(chunk.size() == 1u);

    chunk_proc.on_bytes_written(1);
    BOOST_TEST(chunk_proc.done());
    chunk = chunk_proc.get_chunk(buff);
    BOOST_TEST(chunk.size() == 0u);
}

BOOST_AUTO_TEST_CASE(nonzero_to_size_steps)
{
    chunk_processor chunk_proc;
    chunk_proc.reset(2, 21);  // simulate a previous op
    std::vector<std::uint8_t> buff(10, 0);

    chunk_proc.reset(3, 7);
    BOOST_TEST(!chunk_proc.done());
    auto chunk = chunk_proc.get_chunk(buff);
    BOOST_TEST(chunk.data() == buff.data() + 3);
    BOOST_TEST(chunk.size() == 4u);

    chunk_proc.on_bytes_written(3);
    BOOST_TEST(!chunk_proc.done());
    chunk = chunk_proc.get_chunk(buff);
    BOOST_TEST(chunk.data() == buff.data() + 6);
    BOOST_TEST(chunk.size() == 1u);

    chunk_proc.on_bytes_written(1);
    BOOST_TEST(chunk_proc.done());
    chunk = chunk_proc.get_chunk(buff);
    BOOST_TEST(chunk.size() == 0u);
}

BOOST_AUTO_TEST_CASE(reset_to_empty)
{
    chunk_processor chunk_proc;
    chunk_proc.reset(10, 20);

    chunk_proc.reset();
    BOOST_TEST(chunk_proc.done());
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_CASE(regular_message)
{
    message_writer processor(8);
    std::vector<std::uint8_t> msg_body{0x01, 0x02, 0x03};
    std::uint8_t seqnum = 2;

    // Operation start
    auto mutbuf = processor.prepare_buffer(msg_body.size(), seqnum);
    BOOST_TEST(!processor.done());
    BOOST_TEST(mutbuf.size() == msg_body.size());

    // Simulate serialization
    copy(msg_body, mutbuf);

    // First (and only) chunk
    auto chunk = processor.next_chunk();
    auto expected = create_frame(2, msg_body);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(chunk, expected);

    // On write successful
    processor.on_bytes_written(7);
    BOOST_TEST(seqnum == 3u);
    BOOST_TEST(processor.done());
}

BOOST_AUTO_TEST_CASE(short_writes)
{
    message_writer processor(8);
    std::vector<std::uint8_t> msg_body{0x01, 0x02, 0x03, 0x04, 0x05, 0x06};
    std::uint8_t seqnum = 2;

    // Operation start
    auto mutbuf = processor.prepare_buffer(msg_body.size(), seqnum);
    BOOST_TEST(!processor.done());
    BOOST_TEST(mutbuf.size() == msg_body.size());

    // Simulate serialization
    copy(msg_body, mutbuf);

    // First chunk
    auto chunk = processor.next_chunk();
    auto expected = create_frame(2, msg_body);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(chunk, expected);

    // Write signals partial write
    processor.on_bytes_written(3);
    BOOST_TEST(seqnum == 3u);
    BOOST_TEST(!processor.done());

    // Remaining of the chunk
    chunk = processor.next_chunk();
    span<const std::uint8_t> expected_buff(expected.data() + 3, 7);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(chunk, expected_buff);

    // Another partial write
    processor.on_bytes_written(2);
    BOOST_TEST(seqnum == 3u);
    BOOST_TEST(!processor.done());

    // Remaining of the chunk
    chunk = processor.next_chunk();
    expected_buff = span<const std::uint8_t>(expected.data() + 5, 5);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(chunk, expected_buff);

    // Zero bytes partial writes work correctly
    processor.on_bytes_written(0);
    BOOST_TEST(seqnum == 3u);
    BOOST_TEST(!processor.done());

    // Remaining of the chunk
    chunk = processor.next_chunk();
    expected_buff = span<const std::uint8_t>(expected.data() + 5, 5);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(chunk, expected_buff);

    // Final bytes
    processor.on_bytes_written(5);
    BOOST_TEST(seqnum == 3u);
    BOOST_TEST(processor.done());
}

BOOST_AUTO_TEST_CASE(empty_message)
{
    message_writer processor(8);
    std::uint8_t seqnum = 2;

    // Operation start
    auto mutbuf = processor.prepare_buffer(0, seqnum);
    BOOST_TEST(!processor.done());
    BOOST_TEST(mutbuf.size() == 0u);

    // Chunk should only contain the header
    auto chunk = processor.next_chunk();
    auto expected = create_empty_frame(2);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(chunk, expected);

    // On write successful
    processor.on_bytes_written(4);
    BOOST_TEST(seqnum == 3u);
    BOOST_TEST(processor.done());
}

BOOST_AUTO_TEST_CASE(message_with_max_frame_size_length)
{
    message_writer processor(8);
    std::vector<std::uint8_t> msg_body{0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08};
    std::uint8_t seqnum = 2;

    // Operation start
    auto mutbuf = processor.prepare_buffer(msg_body.size(), seqnum);
    BOOST_TEST(!processor.done());
    BOOST_TEST(mutbuf.size() == msg_body.size());

    // Simulate serialization
    copy(msg_body, mutbuf);

    // Chunk
    auto chunk = processor.next_chunk();
    auto expected = create_frame(2, msg_body);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(chunk, expected);

    // On write successful
    processor.on_bytes_written(12);
    BOOST_TEST(!processor.done());

    // Next chunk is an empty frame
    chunk = processor.next_chunk();
    expected = create_empty_frame(3);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(chunk, expected);

    // On write successful
    processor.on_bytes_written(4);
    BOOST_TEST(seqnum == 4u);
    BOOST_TEST(processor.done());
}

BOOST_AUTO_TEST_CASE(multiframe_message)
{
    message_writer processor(8);
    std::vector<std::uint8_t> msg_frame_1{0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08};
    std::vector<std::uint8_t> msg_frame_2{0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18};
    std::vector<std::uint8_t> msg_frame_3{0x21};
    auto msg = buffer_builder().add(msg_frame_1).add(msg_frame_2).add(msg_frame_3).build();
    std::uint8_t seqnum = 2;

    // Operation start
    auto mutbuf = processor.prepare_buffer(msg.size(), seqnum);
    BOOST_TEST(!processor.done());
    BOOST_TEST(mutbuf.size() == msg.size());

    // Simulate serialization
    copy(msg, mutbuf);

    // Chunk 1
    auto chunk = processor.next_chunk();
    auto expected = create_frame(2, msg_frame_1);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(chunk, expected);

    // On write successful 1 (short write)
    processor.on_bytes_written(4);
    BOOST_TEST(!processor.done());

    // Rest of chunk 1
    chunk = processor.next_chunk();
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(chunk, span<const std::uint8_t>(expected.data() + 4, 8));

    // On write rest of chunk 1
    processor.on_bytes_written(8);
    BOOST_TEST(!processor.done());

    // Chunk 2
    chunk = processor.next_chunk();
    expected = create_frame(3, msg_frame_2);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(chunk, expected);

    // On write successful 2
    processor.on_bytes_written(12);
    BOOST_TEST(!processor.done());

    // Chunk 3
    chunk = processor.next_chunk();
    expected = create_frame(4, msg_frame_3);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(chunk, expected);

    // On write successful
    processor.on_bytes_written(5);
    BOOST_TEST(seqnum == 5u);
    BOOST_TEST(processor.done());
}

BOOST_AUTO_TEST_CASE(multiframe_message_with_max_frame_size)
{
    message_writer processor(8);
    std::vector<std::uint8_t> msg_frame_1{0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08};
    std::vector<std::uint8_t> msg_frame_2{0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18};
    auto msg = concat_copy(msg_frame_1, msg_frame_2);
    std::uint8_t seqnum = 2;

    // Operation start
    auto mutbuf = processor.prepare_buffer(msg.size(), seqnum);
    BOOST_TEST(!processor.done());
    BOOST_TEST(mutbuf.size() == msg.size());

    // Simulate serialization
    copy(msg, mutbuf);

    // Chunk 1
    auto chunk = processor.next_chunk();
    auto expected = create_frame(2, msg_frame_1);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(chunk, expected);

    // On write successful 1
    processor.on_bytes_written(12);
    BOOST_TEST(!processor.done());

    // Chunk 2
    chunk = processor.next_chunk();
    expected = create_frame(3, msg_frame_2);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(chunk, expected);

    // On write successful 2
    processor.on_bytes_written(12);
    BOOST_TEST(!processor.done());

    // Chunk 3 (empty)
    chunk = processor.next_chunk();
    expected = create_empty_frame(4);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(chunk, expected);

    // On write successful
    processor.on_bytes_written(4);
    BOOST_TEST(seqnum == 5u);
    BOOST_TEST(processor.done());
}

BOOST_AUTO_TEST_CASE(seqnum_overflow)
{
    message_writer processor(8);
    std::vector<std::uint8_t> msg{0x01, 0x02};
    std::uint8_t seqnum = 0xff;

    // Operation start
    auto mutbuf = processor.prepare_buffer(msg.size(), seqnum);
    BOOST_TEST(!processor.done());
    BOOST_TEST(mutbuf.size() == msg.size());

    // Simulate serialization
    copy(msg, mutbuf);

    // Prepare chunk
    auto chunk = processor.next_chunk();
    auto expected = create_frame(0xff, msg);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(chunk, expected);

    // On write successful
    processor.on_bytes_written(6);
    BOOST_TEST(seqnum == 0);
    BOOST_TEST(processor.done());
}

BOOST_AUTO_TEST_CASE(several_messages)
{
    message_writer processor(8);
    std::vector<std::uint8_t> msg_1{0x01, 0x02, 0x04};
    std::vector<std::uint8_t> msg_2{0x04, 0x05, 0x06, 0x09, 0xff};
    std::vector<std::uint8_t> msg_3{0x02, 0xab};
    std::uint8_t seqnum_1 = 2;
    std::uint8_t seqnum_2 = 42;
    std::uint8_t seqnum_3 = 21;

    // Operation start for message 1
    auto mutbuf = processor.prepare_buffer(msg_1.size(), seqnum_1);
    BOOST_TEST(!processor.done());
    BOOST_TEST(mutbuf.size() == msg_1.size());
    copy(msg_1, mutbuf);

    // Chunk 1
    auto chunk = processor.next_chunk();
    auto expected = create_frame(2, msg_1);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(chunk, expected);

    // On write successful message 1
    processor.on_bytes_written(7);
    BOOST_TEST(seqnum_1 == 3u);
    BOOST_TEST(processor.done());

    // Operation start for message 2
    mutbuf = processor.prepare_buffer(msg_2.size(), seqnum_2);
    BOOST_TEST(!processor.done());
    BOOST_TEST(mutbuf.size() == msg_2.size());
    copy(msg_2, mutbuf);

    // Chunk 2
    chunk = processor.next_chunk();
    expected = create_frame(42, msg_2);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(chunk, expected);

    // On write successful message 2
    processor.on_bytes_written(9);
    BOOST_TEST(seqnum_2 == 43u);
    BOOST_TEST(processor.done());

    // Operation start for message 3
    mutbuf = processor.prepare_buffer(msg_3.size(), seqnum_3);
    BOOST_TEST(!processor.done());
    BOOST_TEST(mutbuf.size() == msg_3.size());
    copy(msg_3, mutbuf);

    // Chunk 3
    chunk = processor.next_chunk();
    expected = create_frame(21, msg_3);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(chunk, expected);

    // On write successful message 3
    processor.on_bytes_written(6);
    BOOST_TEST(seqnum_3 == 22u);
    BOOST_TEST(processor.done());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace
