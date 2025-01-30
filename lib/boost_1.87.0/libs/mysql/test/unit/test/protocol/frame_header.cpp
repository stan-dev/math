//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/impl/internal/protocol/frame_header.hpp>

#include <boost/core/span.hpp>
#include <boost/test/unit_test.hpp>

#include <cstdint>
#include <memory>

#include "serialization_test.hpp"
#include "test_common/assert_buffer_equals.hpp"

using namespace boost::mysql::detail;
using namespace boost::mysql::test;
using boost::span;

BOOST_AUTO_TEST_CASE(test_frame_header)
{
    struct
    {
        const char* name;
        frame_header header;
        std::array<std::uint8_t, 4> serialized;
    } test_cases[] = {
        {"small_packet_seqnum_0",     {3, 0},           {{0x03, 0x00, 0x00, 0x00}}},
        {"small_packet_seqnum_not_0", {9, 2},           {{0x09, 0x00, 0x00, 0x02}}},
        {"big_packet_seqnum_0",       {0xcacbcc, 0xfa}, {{0xcc, 0xcb, 0xca, 0xfa}}},
        {"max_packet_max_seqnum",     {0xffffff, 0xff}, {{0xff, 0xff, 0xff, 0xff}}}
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name << " serialization")
        {
            std::unique_ptr<std::uint8_t[]> buff{new std::uint8_t[4]};
            span<std::uint8_t, frame_header_size> buff_span(buff.get(), 4);
            serialize_frame_header(buff_span, tc.header);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff_span, tc.serialized);
        }
        BOOST_TEST_CONTEXT(tc.name << " deserialization")
        {
            deserialization_buffer buffer(tc.serialized);
            auto actual = deserialize_frame_header(span<const std::uint8_t, frame_header_size>(buffer));
            BOOST_TEST(actual.size == tc.header.size);
            BOOST_TEST(actual.sequence_number == tc.header.sequence_number);
        }
    }
}
