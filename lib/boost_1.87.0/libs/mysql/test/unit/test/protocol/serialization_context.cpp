//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/impl/internal/protocol/impl/serialization_context.hpp>

#include <boost/test/unit_test.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/buffer_concat.hpp"
#include "test_common/printing.hpp"

using namespace boost::mysql;

namespace {

BOOST_AUTO_TEST_SUITE(test_serialization_context)

struct framing_test_case
{
    string_view name;
    std::size_t expected_next_frame_offset;  // not counting previous contents
    std::vector<std::uint8_t> payload;
    std::vector<std::uint8_t> expected_buffer;  // not counting previous contents
};

std::vector<framing_test_case> make_test_cases()
{
    return {
        // clang-format off
        {"0 bytes",     12,   {},                          {0, 0, 0, 0}                                       },
        {"1 byte",      12,   {1},                         {0, 0, 0, 0, 1}                                    },
        {"5 bytes",     12,   {1, 2, 3, 4, 5},             {0, 0, 0, 0, 1, 2, 3, 4, 5}                        },
        {"fs-1 bytes",  12,   {1, 2, 3, 4, 5, 6, 7},       {0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7}                  },
        {"fs bytes",    24,   {1, 2, 3, 4, 5, 6, 7, 8},    {0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0}   },
        {"fs+1 bytes",  24,   {1, 2, 3, 4, 5, 6, 7, 8, 9}, {0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 9}},
        {"2fs bytes",   36,   {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16},
            {0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 9, 10, 11, 12, 13, 14, 15, 16, 0, 0, 0, 0}    },
        {"2fs+1 bytes", 36,   {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17},
            {0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 9, 10, 11, 12, 13, 14, 15, 16, 0, 0, 0, 0, 17}},
        {"2fs+5 bytes", 36,   {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21},
            {0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 9, 10, 11, 12, 13, 14, 15, 16, 0, 0, 0, 0, 17, 18, 19, 20, 21}},
        {"3fs-1 bytes", 36,   {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23},
            {0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 9, 10, 11, 12, 13, 14, 15, 16, 0, 0, 0, 0, 17, 18, 19, 20, 21, 22, 23}},
        {"3fs bytes",   48,   {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24},
            {0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 9, 10, 11, 12, 13, 14, 15, 16, 0, 0, 0, 0, 17, 18, 19, 20, 21, 22, 23, 24, 0, 0, 0, 0}},
        {"3fs+1 bytes", 48,   {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25},
            {0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 9, 10, 11, 12, 13, 14, 15, 16, 0, 0, 0, 0, 17, 18, 19, 20, 21, 22, 23, 24, 0, 0, 0, 0, 25}},
        // clang-format on
    };
}

BOOST_AUTO_TEST_CASE(add)
{
    constexpr std::size_t fs = 8u;  // frame size
    const std::vector<std::uint8_t> initial_buffer{0xaa, 0xbb, 0xcc, 0xdd, 0xee};

    for (const auto& tc : make_test_cases())
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // Setup
            std::vector<std::uint8_t> buff{initial_buffer};
            detail::serialization_context ctx(buff, 0xffff, fs);

            // Add the payload
            ctx.add(tc.payload);

            // Check
            auto expected = test::concat(initial_buffer, tc.expected_buffer);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);
            BOOST_TEST(ctx.next_header_offset() == tc.expected_next_frame_offset + initial_buffer.size());
            BOOST_TEST(ctx.error() == error_code());
        }
    }
}

// Spotcheck: if the initial buffer is empty, everything works fine
BOOST_AUTO_TEST_CASE(add_initial_buffer_empty)
{
    // Setup
    std::vector<std::uint8_t> buff;
    detail::serialization_context ctx(buff, 0xffff, 8);

    // Add data
    const std::array<std::uint8_t, 10> payload{
        {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
    };
    ctx.add(payload);

    // Check
    const std::vector<std::uint8_t> expected{0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 9, 10};
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);
    BOOST_TEST(ctx.next_header_offset() == 24u);
    BOOST_TEST(ctx.error() == error_code());
}

// Spotcheck: adding single bytes or in chunks also works fine
BOOST_AUTO_TEST_CASE(chunks)
{
    // Setup
    std::vector<std::uint8_t> buff;
    detail::serialization_context ctx(buff, 0xffff, 8);
    const std::array<std::uint8_t, 4> payload1{
        {1, 2, 3, 4}
    };
    const std::array<std::uint8_t, 5> payload2{
        {5, 6, 7, 8, 9}
    };

    // Add byte
    ctx.add(0xff);
    std::vector<std::uint8_t> expected{0, 0, 0, 0, 0xff};
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);
    BOOST_TEST(ctx.error() == error_code());

    // Add buffer
    ctx.add(payload1);
    expected = {0, 0, 0, 0, 0xff, 1, 2, 3, 4};
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);
    BOOST_TEST(ctx.error() == error_code());

    // Add byte
    ctx.add(0xfe);
    expected = {0, 0, 0, 0, 0xff, 1, 2, 3, 4, 0xfe};
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);
    BOOST_TEST(ctx.error() == error_code());

    // Add buffer
    ctx.add(payload2);
    expected = {0, 0, 0, 0, 0xff, 1, 2, 3, 4, 0xfe, 5, 6, 0, 0, 0, 0, 7, 8, 9};
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);
    BOOST_TEST(ctx.error() == error_code());

    // Add byte
    ctx.add(0xfc);
    expected = {0, 0, 0, 0, 0xff, 1, 2, 3, 4, 0xfe, 5, 6, 0, 0, 0, 0, 7, 8, 9, 0xfc};
    BOOST_TEST(ctx.next_header_offset() == 24u);
    BOOST_TEST(ctx.error() == error_code());
}

// Spotcheck: adding a single byte that causes a frame header to be written works
BOOST_AUTO_TEST_CASE(add_byte_fills_frame)
{
    // Setup
    std::vector<std::uint8_t> buff;
    detail::serialization_context ctx(buff, 0xffff, 8);
    const std::array<std::uint8_t, 7> payload{
        {1, 2, 3, 4, 5, 6, 7}
    };

    // Add payload
    ctx.add(payload);
    std::vector<std::uint8_t> expected{0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7};
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);
    BOOST_TEST(ctx.next_header_offset() == 12u);
    BOOST_TEST(ctx.error() == error_code());

    // Add byte
    ctx.add(0xab);
    expected = {0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 0xab, 0, 0, 0, 0};
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);
    BOOST_TEST(ctx.next_header_offset() == 24u);
    BOOST_TEST(ctx.error() == error_code());
}

BOOST_AUTO_TEST_CASE(write_frame_headers)
{
    struct
    {
        string_view name;
        std::uint8_t expected_seqnum;
        std::vector<std::uint8_t> payload;
        std::vector<std::uint8_t> expected;
    } test_cases[] = {
        // clang-format off
        {"0 bytes",     43,   {},                          {0, 0, 0, 42}                                       },
        {"1 byte",      43,   {1},                         {1, 0, 0, 42, 1}                                    },
        {"5 bytes",     43,   {1, 2, 3, 4, 5},             {5, 0, 0, 42, 1, 2, 3, 4, 5}                        },
        {"fs-1 bytes",  43,   {1, 2, 3, 4, 5, 6, 7},       {7, 0, 0, 42, 1, 2, 3, 4, 5, 6, 7}                  },
        {"fs bytes",    44,   {1, 2, 3, 4, 5, 6, 7, 8},    {8, 0, 0, 42, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 43}   },
        {"fs+1 bytes",  44,   {1, 2, 3, 4, 5, 6, 7, 8, 9}, {8, 0, 0, 42, 1, 2, 3, 4, 5, 6, 7, 8, 1, 0, 0, 43, 9}},
        {"2fs bytes",   45,   {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16},
            {8, 0, 0, 42, 1, 2, 3, 4, 5, 6, 7, 8, 8, 0, 0, 43, 9, 10, 11, 12, 13, 14, 15, 16, 0, 0, 0, 44}    },
        {"2fs+1 bytes", 45,   {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17},
            {8, 0, 0, 42, 1, 2, 3, 4, 5, 6, 7, 8, 8, 0, 0, 43, 9, 10, 11, 12, 13, 14, 15, 16, 1, 0, 0, 44, 17}},
        {"2fs+5 bytes", 45,   {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21},
            {8, 0, 0, 42, 1, 2, 3, 4, 5, 6, 7, 8, 8, 0, 0, 43, 9, 10, 11, 12, 13, 14, 15, 16, 5, 0, 0, 44, 17, 18, 19, 20, 21}},
        {"3fs-1 bytes", 45,   {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23},
            {8, 0, 0, 42, 1, 2, 3, 4, 5, 6, 7, 8, 8, 0, 0, 43, 9, 10, 11, 12, 13, 14, 15, 16, 7, 0, 0, 44, 17, 18, 19, 20, 21, 22, 23}},
        {"3fs bytes",   46,   {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24},
            {8, 0, 0, 42, 1, 2, 3, 4, 5, 6, 7, 8, 8, 0, 0, 43, 9, 10, 11, 12, 13, 14, 15, 16, 8, 0, 0, 44, 17, 18, 19, 20, 21, 22, 23, 24,
                0, 0, 0, 45}},
        {"3fs+1 bytes", 46,   {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25},
            {8, 0, 0, 42, 1, 2, 3, 4, 5, 6, 7, 8, 8, 0, 0, 43, 9, 10, 11, 12, 13, 14, 15, 16, 8, 0, 0, 44, 17, 18, 19, 20, 21, 22, 23, 24,
                1, 0, 0, 45, 25}},
        // clang-format on
    };

    const std::vector<std::uint8_t> initial_buffer{90, 91, 92, 93, 94};

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // Setup
            std::vector<std::uint8_t> buff{initial_buffer};
            detail::serialization_context ctx(buff, 0xffff, 8);
            ctx.add(tc.payload);

            // Call and check
            auto seqnum = ctx.write_frame_headers(42, initial_buffer.size());
            const auto expected = test::concat(initial_buffer, tc.expected);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);
            BOOST_TEST(seqnum == tc.expected_seqnum);
            BOOST_TEST(ctx.error() == error_code());
        }
    }
}

// Spotcheck: we correctly wrap sequence numbers when going over 0xff
BOOST_AUTO_TEST_CASE(write_frame_headers_seqnum_wrap)
{
    // Setup
    std::vector<std::uint8_t> buff;
    detail::serialization_context ctx(buff, 0xffff, 8);
    for (std::uint8_t i = 1; i <= 20; ++i)
        ctx.add(i);

    // Call and check
    const std::vector<std::uint8_t> expected{
        8, 0, 0, 0xfe, 1,  2,  3,  4,  5,  6,  7,  8,   // frame 1
        8, 0, 0, 0xff, 9,  10, 11, 12, 13, 14, 15, 16,  // frame 2
        4, 0, 0, 0,    17, 18, 19, 20                   // frame 3
    };
    auto seqnum = ctx.write_frame_headers(0xfe, 0);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);
    BOOST_TEST(seqnum == 1u);
    BOOST_TEST(ctx.error() == error_code());
}

// Spotcheck: disable framing works
BOOST_AUTO_TEST_CASE(disable_framing)
{
    // Setup
    std::vector<std::uint8_t> buff;
    detail::serialization_context ctx(buff, 0xffff, detail::disable_framing);

    // Add data using the several functions available
    const std::array<std::uint8_t, 5> payload1{
        {1, 2, 3, 4, 5}
    };
    const std::array<std::uint8_t, 4> payload2{
        {6, 7, 8, 9}
    };
    ctx.add(42);
    ctx.add(payload1);
    ctx.add(payload2);

    // We didn't add any framing
    const std::vector<std::uint8_t> expected{42, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);
    BOOST_TEST(ctx.error() == error_code());
}

BOOST_AUTO_TEST_SUITE(max_buffer_size_error)

BOOST_AUTO_TEST_CASE(header_exceeds_maxsize)
{
    // Setup
    std::vector<std::uint8_t> buff;
    detail::serialization_context ctx(buff, 3u);

    // Buffer can't hold the header
    BOOST_TEST(ctx.error() == client_errc::max_buffer_size_exceeded);
}

BOOST_AUTO_TEST_CASE(contents_exceed_maxsize)
{
    // Setup
    std::vector<std::uint8_t> buff;
    detail::serialization_context ctx(buff, 8u);

    // Header plus content would exceed max size
    ctx.add(std::vector<std::uint8_t>{1, 2, 3, 4, 5, 6});
    BOOST_TEST(ctx.error() == client_errc::max_buffer_size_exceeded);

    // Only header written
    std::array<std::uint8_t, 4> expected{};
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);
}

BOOST_AUTO_TEST_CASE(subsequent_header_exceeds_maxsize)
{
    // Setup
    std::vector<std::uint8_t> buff;
    detail::serialization_context ctx(buff, 13u, 8u);

    // Successfully add some data
    ctx.add(std::vector<std::uint8_t>{1, 2, 3, 4, 5, 6});
    std::vector<std::uint8_t> expected{0, 0, 0, 0, 1, 2, 3, 4, 5, 6};
    BOOST_TEST(ctx.error() == error_code());
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);

    // Add data triggering a header that can't fit
    ctx.add(std::vector<std::uint8_t>{7, 8});
    BOOST_TEST(ctx.error() == client_errc::max_buffer_size_exceeded);
}

BOOST_AUTO_TEST_CASE(maxsize_zero)
{
    // Setup
    std::vector<std::uint8_t> buff;
    detail::serialization_context ctx(buff, 0u);

    // Buffer can't hold the header
    BOOST_TEST(ctx.error() == client_errc::max_buffer_size_exceeded);
}

BOOST_AUTO_TEST_CASE(error_by_one_byte)
{
    // Setup
    std::vector<std::uint8_t> buff;
    detail::serialization_context ctx(buff, 12u);

    // Successfully add data until max size
    ctx.add(std::vector<std::uint8_t>{1, 2, 3, 4, 5, 6, 7, 8});
    std::vector<std::uint8_t> expected{0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8};
    BOOST_TEST(ctx.error() == error_code());
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);

    // Adding more data fails. No data is written to the buffer
    ctx.add(std::vector<std::uint8_t>{1});
    BOOST_TEST(ctx.error() == client_errc::max_buffer_size_exceeded);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);
}

// Spotcheck: adding a single byte triggers the same behavior
BOOST_AUTO_TEST_CASE(error_add_u8)
{
    // Setup
    std::vector<std::uint8_t> buff;
    detail::serialization_context ctx(buff, 12u);

    // Successfully add data until max size
    ctx.add(std::vector<std::uint8_t>{1, 2, 3, 4, 5, 6, 7, 8});
    std::vector<std::uint8_t> expected{0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8};
    BOOST_TEST(ctx.error() == error_code());
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);

    // Adding more data fails. No data is written to the buffer
    ctx.add(42);
    BOOST_TEST(ctx.error() == client_errc::max_buffer_size_exceeded);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);
}

// Edge case: if the input buffer exceeded the max size, we fail
BOOST_AUTO_TEST_CASE(buffer_exceeds_max_size)
{
    std::vector<std::uint8_t> buff(48u, 0);
    detail::serialization_context ctx(buff, 12u);
    BOOST_TEST(ctx.error() == client_errc::max_buffer_size_exceeded);
}

// Previous contents are taken into account for size checks
BOOST_AUTO_TEST_CASE(buffer_with_previous_contents)
{
    // Setup
    std::vector<std::uint8_t> buff{1, 2, 3};
    detail::serialization_context ctx(buff, 8u);

    // Just max size
    ctx.add(42);
    std::vector<std::uint8_t> expected{1, 2, 3, 0, 0, 0, 0, 42};
    BOOST_TEST(ctx.error() == error_code());
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);

    // Past max size
    ctx.add(80);
    BOOST_TEST(ctx.error() == client_errc::max_buffer_size_exceeded);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);
}

BOOST_AUTO_TEST_CASE(several_errors)
{
    // Setup
    std::vector<std::uint8_t> buff;
    detail::serialization_context ctx(buff, 12u);

    // Successfully add some data
    ctx.add(std::vector<std::uint8_t>{1, 2, 3, 4, 5});
    std::vector<std::uint8_t> expected{0, 0, 0, 0, 1, 2, 3, 4, 5};
    BOOST_TEST(ctx.error() == error_code());
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);

    // Adding more data fails. No data is written to the buffer
    ctx.add(std::vector<std::uint8_t>{6, 7, 8, 9});
    BOOST_TEST(ctx.error() == client_errc::max_buffer_size_exceeded);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);

    // Adding more data again does nothing
    ctx.add(std::vector<std::uint8_t>{10, 11, 12});
    BOOST_TEST(ctx.error() == client_errc::max_buffer_size_exceeded);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);
}

BOOST_AUTO_TEST_CASE(success_after_error)
{
    // Setup
    std::vector<std::uint8_t> buff;
    detail::serialization_context ctx(buff, 12u);

    // Successfully add some data
    ctx.add(std::vector<std::uint8_t>{1, 2, 3, 4, 5});
    std::vector<std::uint8_t> expected{0, 0, 0, 0, 1, 2, 3, 4, 5};
    BOOST_TEST(ctx.error() == error_code());
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);

    // Adding more data fails. No data is written to the buffer
    ctx.add(std::vector<std::uint8_t>{6, 7, 8, 9});
    BOOST_TEST(ctx.error() == client_errc::max_buffer_size_exceeded);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);

    // Adding more data again does nothing, even if the data would fit
    ctx.add(1);
    BOOST_TEST(ctx.error() == client_errc::max_buffer_size_exceeded);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);
}

BOOST_AUTO_TEST_CASE(add_error)
{
    // Setup
    std::vector<std::uint8_t> buff;
    detail::serialization_context ctx(buff);

    // Successfully add some data
    ctx.add(std::vector<std::uint8_t>{1, 2, 3, 4, 5});
    std::vector<std::uint8_t> expected{0, 0, 0, 0, 1, 2, 3, 4, 5};
    BOOST_TEST(ctx.error() == error_code());
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);

    // Add an error
    ctx.add_error(client_errc::invalid_encoding);
    BOOST_TEST(ctx.error() == client_errc::invalid_encoding);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);

    // Adding further data with the error set does nothing
    ctx.add(std::vector<std::uint8_t>{8, 9, 10});
    BOOST_TEST(ctx.error() == client_errc::invalid_encoding);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);

    // Adding another error does nothing
    ctx.add_error(client_errc::protocol_value_error);
    BOOST_TEST(ctx.error() == client_errc::invalid_encoding);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);

    ctx.add_error(error_code());
    BOOST_TEST(ctx.error() == client_errc::invalid_encoding);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

}  // namespace