//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/impl/internal/channel/message_parser.hpp>
#include <boost/mysql/impl/internal/channel/read_buffer.hpp>

#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/buffer_concat.hpp"
#include "test_unit/create_frame.hpp"

using namespace boost::mysql::detail;
using namespace boost::mysql::test;
using boost::span;

BOOST_AUTO_TEST_SUITE(test_message_parser)

class parser_fixture
{
    message_parser parser_{64};  // max frame size is 64
    read_buffer buff_;
    std::vector<std::uint8_t> contents_;
    std::size_t bytes_written_{0};
    const std::uint8_t* buffer_first_;

public:
    parser_fixture(std::vector<std::uint8_t> contents, std::size_t buffsize = 512)
        : buff_(buffsize), contents_(std::move(contents)), buffer_first_(buff_.first())
    {
    }
    read_buffer& buffer() noexcept { return buff_; }
    message_parser::result parse_bytes(std::size_t num_bytes)
    {
        message_parser::result res;
        if (num_bytes)
        {
            std::memcpy(buff_.free_first(), contents_.data() + bytes_written_, num_bytes);
            bytes_written_ += num_bytes;
            buff_.move_to_pending(num_bytes);
        }
        parser_.parse_message(buff_, res);
        return res;
    }
    span<const std::uint8_t> check_message(const std::vector<std::uint8_t>& contents)
    {
        span<const std::uint8_t> msg(buff_.current_message_first() - contents.size(), contents.size());
        BOOST_MYSQL_ASSERT_BUFFER_EQUALS(msg, contents);
        return msg;
    }
    void check_buffer_stability() { BOOST_TEST(buff_.first() == buffer_first_); }
};

BOOST_AUTO_TEST_CASE(fragmented_header_and_body_multiple)
{
    // message to be parsed
    parser_fixture fixture(create_frame(0, {0x01, 0x02, 0x03}));

    // 1 byte in the header received
    auto res = fixture.parse_bytes(1);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 3u);

    // Another 2 bytes received
    res = fixture.parse_bytes(2);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 1u);

    // Header fully received
    res = fixture.parse_bytes(1);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 3u);

    // 1 byte in body received
    res = fixture.parse_bytes(1);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 2u);

    // body fully received (single frame messages keep header as an optimization)
    res = fixture.parse_bytes(2);
    fixture.check_message({0x01, 0x02, 0x03});
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 3u);
    BOOST_TEST(res.message.seqnum_first == 0u);
    BOOST_TEST(res.message.seqnum_last == 0u);
    BOOST_TEST(!res.message.has_seqnum_mismatch);

    // Buffer did not reallocate
    fixture.check_buffer_stability();
}

BOOST_AUTO_TEST_CASE(fragmented_header_and_body_single)
{
    // message to be parsed
    parser_fixture fixture(create_frame(0, {0x01, 0x02, 0x03}));

    // Full header received
    auto res = fixture.parse_bytes(4);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 3u);

    // Full body received
    res = fixture.parse_bytes(3);
    fixture.check_message({0x01, 0x02, 0x03});
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 3u);
    BOOST_TEST(res.message.seqnum_first == 0u);
    BOOST_TEST(res.message.seqnum_last == 0u);
    BOOST_TEST(!res.message.has_seqnum_mismatch);

    // Buffer did not reallocate
    fixture.check_buffer_stability();
}

BOOST_AUTO_TEST_CASE(fragmented_body)
{
    // message to be parsed
    parser_fixture fixture(create_frame(0, {0x01, 0x02, 0x03}));

    // Full header and body part received
    auto res = fixture.parse_bytes(5);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 2u);

    // Full body + next header + next body received
    res = fixture.parse_bytes(2);
    fixture.check_message({0x01, 0x02, 0x03});
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 3u);
    BOOST_TEST(res.message.seqnum_first == 0u);
    BOOST_TEST(res.message.seqnum_last == 0u);
    BOOST_TEST(!res.message.has_seqnum_mismatch);

    // Buffer did not reallocate
    fixture.check_buffer_stability();
}

BOOST_AUTO_TEST_CASE(full_message)
{
    // message to be parsed
    parser_fixture fixture(create_frame(0, {0x01, 0x02, 0x03}));

    // Full header and body part received
    auto res = fixture.parse_bytes(7);
    fixture.check_message({0x01, 0x02, 0x03});
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 3u);
    BOOST_TEST(res.message.seqnum_first == 0u);
    BOOST_TEST(res.message.seqnum_last == 0u);
    BOOST_TEST(!res.message.has_seqnum_mismatch);

    // Buffer did not reallocate
    fixture.check_buffer_stability();
}

BOOST_AUTO_TEST_CASE(empty_message)
{
    // message to be parsed
    parser_fixture fixture(create_empty_frame(1));

    // Full header and body part received
    auto res = fixture.parse_bytes(4);
    fixture.check_message({});
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 0u);
    BOOST_TEST(res.message.seqnum_first == 1u);
    BOOST_TEST(res.message.seqnum_last == 1u);
    BOOST_TEST(!res.message.has_seqnum_mismatch);

    // Buffer did not reallocate
    fixture.check_buffer_stability();
}

BOOST_AUTO_TEST_CASE(two_messages_one_after_another)
{
    // message to be parsed
    std::vector<std::uint8_t> first_msg_body{0x01, 0x02, 0x03};
    std::vector<std::uint8_t> second_msg_body{0x04, 0x05, 0x06, 0x07};
    parser_fixture fixture(concat_copy(create_frame(0, first_msg_body), create_frame(2, second_msg_body)));

    // 1st message
    auto res = fixture.parse_bytes(7);
    auto first_msg = fixture.check_message(first_msg_body);
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 3u);
    BOOST_TEST(res.message.seqnum_first == 0u);
    BOOST_TEST(res.message.seqnum_last == 0u);
    BOOST_TEST(!res.message.has_seqnum_mismatch);

    // 2nd message
    res = fixture.parse_bytes(8);
    fixture.check_message(second_msg_body);
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 4u);
    BOOST_TEST(res.message.seqnum_first == 2u);
    BOOST_TEST(res.message.seqnum_last == 2u);
    BOOST_TEST(!res.message.has_seqnum_mismatch);

    // 1st message still valid
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(first_msg, first_msg_body);

    // Buffer did not reallocate
    fixture.check_buffer_stability();
}

BOOST_AUTO_TEST_CASE(two_messages_at_once)
{
    // message to be parsed
    std::vector<std::uint8_t> first_msg_body{0x01, 0x02, 0x03};
    std::vector<std::uint8_t> second_msg_body{0x04, 0x05, 0x06, 0x07};
    parser_fixture fixture(concat_copy(create_frame(0, first_msg_body), create_frame(2, second_msg_body)));

    // 1st message
    auto res = fixture.parse_bytes(15);
    auto first_msg = fixture.check_message(first_msg_body);
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 3u);
    BOOST_TEST(res.message.seqnum_first == 0u);
    BOOST_TEST(res.message.seqnum_last == 0u);
    BOOST_TEST(!res.message.has_seqnum_mismatch);

    // 2nd message
    res = fixture.parse_bytes(0);
    fixture.check_message(second_msg_body);
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 4u);
    BOOST_TEST(res.message.seqnum_first == 2u);
    BOOST_TEST(res.message.seqnum_last == 2u);
    BOOST_TEST(!res.message.has_seqnum_mismatch);

    // 1st message still valid
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(first_msg, first_msg_body);

    // Buffer did not reallocate
    fixture.check_buffer_stability();
}

BOOST_AUTO_TEST_CASE(three_messages_last_fragmented)
{
    // message to be parsed
    std::vector<std::uint8_t> first_msg_body{0x01, 0x02, 0x03};
    std::vector<std::uint8_t> second_msg_body{0x04, 0x05, 0x06, 0x07};
    std::vector<std::uint8_t> third_msg_body{0x08, 0x09};
    parser_fixture fixture(buffer_builder()
                               .add(create_frame(0, first_msg_body))
                               .add(create_frame(2, second_msg_body))
                               .add(create_frame(3, third_msg_body))
                               .build());

    // 1st message
    auto res = fixture.parse_bytes(20);  // 1st and 2nd messages + 3rd message header and body part
    auto first_msg = fixture.check_message(first_msg_body);
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 3u);
    BOOST_TEST(res.message.seqnum_first == 0u);
    BOOST_TEST(res.message.seqnum_last == 0u);
    BOOST_TEST(!res.message.has_seqnum_mismatch);

    // 2nd message
    res = fixture.parse_bytes(0);
    auto second_msg = fixture.check_message(second_msg_body);
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 4u);
    BOOST_TEST(res.message.seqnum_first == 2u);
    BOOST_TEST(res.message.seqnum_last == 2u);
    BOOST_TEST(!res.message.has_seqnum_mismatch);

    // 3rd message header + body part
    res = fixture.parse_bytes(0);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 1u);

    // 3rd message
    res = fixture.parse_bytes(1);
    fixture.check_message(third_msg_body);
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 2u);
    BOOST_TEST(res.message.seqnum_first == 3u);
    BOOST_TEST(res.message.seqnum_last == 3u);
    BOOST_TEST(!res.message.has_seqnum_mismatch);

    // 1st and 2nd messages still valid
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(first_msg, first_msg_body);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(second_msg, second_msg_body);

    // Buffer did not reallocate
    fixture.check_buffer_stability();
}

BOOST_AUTO_TEST_CASE(two_frame_message)
{
    // message to be parsed
    parser_fixture fixture(
        concat_copy(
            create_frame(0, std::vector<std::uint8_t>(64, 0x04)),
            create_frame(1, {0x05, 0x06, 0x07})
        ),
        64 + 16
    );
    auto expected_message = concat_copy(std::vector<std::uint8_t>(64, 0x04), {0x05, 0x06, 0x07});

    // header 1 + body part
    auto res = fixture.parse_bytes(6);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == unsigned(64 - 2));

    // body part
    res = fixture.parse_bytes(64 - 10);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 8u);

    // body part + header 2 + body 2 part
    res = fixture.parse_bytes(13);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 2u);

    // remaining of body 2
    res = fixture.parse_bytes(2);
    fixture.check_message(expected_message);
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 64u + 3u);
    BOOST_TEST(res.message.seqnum_first == 0u);
    BOOST_TEST(res.message.seqnum_last == 1u);
    BOOST_TEST(!res.message.has_seqnum_mismatch);

    // buffer did not reallocate
    fixture.check_buffer_stability();
}

BOOST_AUTO_TEST_CASE(two_frame_message_with_reserved_area)
{
    // message to be parsed
    std::vector<std::uint8_t> first_msg_body{0x01, 0x02, 0x03};
    parser_fixture fixture(
        buffer_builder()
            .add(create_frame(0, first_msg_body))
            .add(create_frame(4, std::vector<std::uint8_t>(64, 0x04)))
            .add(create_frame(5, {0x05, 0x06, 0x07}))
            .build(),
        64 + 64
    );
    auto second_msg_body = concat_copy(std::vector<std::uint8_t>(64, 0x04), {0x05, 0x06, 0x07});

    // msg 1
    auto res = fixture.parse_bytes(7);
    auto first_msg = fixture.check_message(first_msg_body);
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 3u);
    BOOST_TEST(res.message.seqnum_first == 0u);
    BOOST_TEST(res.message.seqnum_last == 0u);
    BOOST_TEST(!res.message.has_seqnum_mismatch);

    // msg 2 (multiframe)
    res = fixture.parse_bytes(64 + 4 * 2 + 3);
    fixture.check_message(second_msg_body);
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 64u + 3u);
    BOOST_TEST(res.message.seqnum_first == 4u);
    BOOST_TEST(res.message.seqnum_last == 5u);
    BOOST_TEST(!res.message.has_seqnum_mismatch);

    // msg 1 still valid
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(first_msg, first_msg_body);

    // buffer did not reallocate
    fixture.check_buffer_stability();
}

BOOST_AUTO_TEST_CASE(two_frame_message_fragmented)
{
    // message to be parsed
    parser_fixture fixture(
        buffer_builder()
            .add(create_frame(0, std::vector<std::uint8_t>(64, 0x04)))
            .add(create_frame(1, {0x05, 0x06, 0x07}))
            .build(),
        64 + 16
    );
    auto expected_message = concat_copy(std::vector<std::uint8_t>(64, 0x04), {0x05, 0x06, 0x07});

    // part of header 1
    auto res = fixture.parse_bytes(3);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 1u);

    // header 1 full
    res = fixture.parse_bytes(1);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 64u);

    // part of body 1
    res = fixture.parse_bytes(64 - 8);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 8u);

    // rest of body 1
    res = fixture.parse_bytes(8);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 4u);

    // part of header 2
    res = fixture.parse_bytes(1);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 3u);

    // another part of header 2
    res = fixture.parse_bytes(2);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 1u);

    // rest of header 2
    res = fixture.parse_bytes(1);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 3u);

    // part of body 1
    res = fixture.parse_bytes(1);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 2u);

    // another part of body 2
    res = fixture.parse_bytes(1);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 1u);

    // remaining of body 2
    res = fixture.parse_bytes(1);
    fixture.check_message(expected_message);
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 64u + 3u);
    BOOST_TEST(res.message.seqnum_first == 0u);
    BOOST_TEST(res.message.seqnum_last == 1u);
    BOOST_TEST(!res.message.has_seqnum_mismatch);

    // buffer did not reallocate
    fixture.check_buffer_stability();
}

BOOST_AUTO_TEST_CASE(three_frame_message)
{
    // message to be parsed
    parser_fixture fixture(
        buffer_builder()
            .add(create_frame(2, std::vector<std::uint8_t>(64, 0x04)))
            .add(create_frame(3, std::vector<std::uint8_t>(64, 0x05)))
            .add(create_frame(4, {0x05, 0x06, 0x07}))
            .build(),
        64 * 2 + 64
    );
    auto expected_message = buffer_builder()
                                .add(std::vector<std::uint8_t>(64, 0x04))
                                .add(std::vector<std::uint8_t>(64, 0x05))
                                .add({0x05, 0x06, 0x07})
                                .build();

    // header 1 + body 1 part
    auto res = fixture.parse_bytes(6);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == unsigned(64u - 2u));

    // body 1 part + header 2 + body 2 part
    res = fixture.parse_bytes(64 - 2 + 4 + 3);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 64u - 3u);

    // body 2 part + header 3 + body 3 part
    res = fixture.parse_bytes(64 - 3 + 4 + 1);
    BOOST_TEST(!res.has_message);
    BOOST_TEST(res.required_size == 2u);

    // body 3 part
    res = fixture.parse_bytes(2);
    fixture.check_message(expected_message);
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 64u * 2u + 3u);
    BOOST_TEST(res.message.seqnum_first == 2u);
    BOOST_TEST(res.message.seqnum_last == 4u);
    BOOST_TEST(!res.message.has_seqnum_mismatch);

    // buffer did not reallocate
    fixture.check_buffer_stability();
}

BOOST_AUTO_TEST_CASE(two_frame_message_mismatched_seqnums)
{
    // message to be parsed
    parser_fixture fixture(
        buffer_builder()
            .add(create_frame(1, std::vector<std::uint8_t>(64, 0x04)))
            .add(create_frame(3, {0x05, 0x06, 0x07}))
            .build(),
        64 + 16
    );
    auto expected_message = concat_copy(std::vector<std::uint8_t>(64, 0x04), {0x05, 0x06, 0x07});

    // all in one
    auto res = fixture.parse_bytes(64 + 4 * 2 + 3);
    fixture.check_message(expected_message);
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 64u + 3u);
    BOOST_TEST(res.message.seqnum_first == 1u);
    BOOST_TEST(res.message.has_seqnum_mismatch);

    // buffer did not reallocate
    fixture.check_buffer_stability();
}

BOOST_AUTO_TEST_CASE(three_frame_message_mismatched_seqnums)
{
    // message to be parsed
    parser_fixture fixture(
        buffer_builder()
            .add(create_frame(1, std::vector<std::uint8_t>(64, 0x04)))
            .add(create_frame(2, std::vector<std::uint8_t>(64, 0x05)))
            .add(create_frame(0, {0x05, 0x06, 0x07}))
            .build(),
        64 * 2 + 64
    );
    auto expected_message = buffer_builder()
                                .add(std::vector<std::uint8_t>(64, 0x04))
                                .add(std::vector<std::uint8_t>(64, 0x05))
                                .add({0x05, 0x06, 0x07})
                                .build();

    // all in one
    auto res = fixture.parse_bytes(64 * 2 + 4 * 3 + 3);
    fixture.check_message(expected_message);
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 64u * 2u + 3u);
    BOOST_TEST(res.message.seqnum_first == 1u);
    BOOST_TEST(res.message.has_seqnum_mismatch);

    // buffer did not reallocate
    fixture.check_buffer_stability();
}

BOOST_AUTO_TEST_CASE(two_frame_seqnum_overflow)
{
    // message to be parsed
    parser_fixture fixture(
        buffer_builder()
            .add(create_frame(255, std::vector<std::uint8_t>(64, 0x04)))
            .add(create_frame(0, {0x05, 0x06, 0x07}))
            .build(),
        64 + 16
    );
    auto expected_message = concat_copy(std::vector<std::uint8_t>(64, 0x04), {0x05, 0x06, 0x07});

    // all in one
    auto res = fixture.parse_bytes(64 + 4 * 2 + 3);
    fixture.check_message(expected_message);
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 64u + 3u);
    BOOST_TEST(res.message.seqnum_first == 255u);
    BOOST_TEST(res.message.seqnum_last == 0u);
    BOOST_TEST(!res.message.has_seqnum_mismatch);

    // buffer did not reallocate
    fixture.check_buffer_stability();
}

BOOST_AUTO_TEST_CASE(two_frame_max_size)
{
    // message to be parsed. The two frames have size == max_frame_size,
    // so a third, empty header is received
    parser_fixture fixture(
        buffer_builder()
            .add(create_frame(1, std::vector<std::uint8_t>(64, 0x04)))
            .add(create_frame(2, std::vector<std::uint8_t>(64, 0x05)))
            .add(create_empty_frame(3))
            .build(),
        64 * 3
    );
    auto expected_message = concat_copy(
        std::vector<std::uint8_t>(64, 0x04),
        std::vector<std::uint8_t>(64, 0x05)
    );

    // all in one
    auto res = fixture.parse_bytes(64 * 2 + 4 * 3);
    fixture.check_message(expected_message);
    BOOST_TEST(res.has_message);
    BOOST_TEST(res.message.size == 64u * 2u);
    BOOST_TEST(res.message.seqnum_first == 1u);
    BOOST_TEST(res.message.seqnum_last == 3u);
    BOOST_TEST(!res.message.has_seqnum_mismatch);

    // buffer did not reallocate
    fixture.check_buffer_stability();
}

BOOST_AUTO_TEST_SUITE_END()
