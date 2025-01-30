//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/field_view.hpp>

#include <boost/mysql/impl/internal/protocol/impl/null_bitmap.hpp>

#include <boost/core/span.hpp>
#include <boost/test/unit_test.hpp>

#include <array>
#include <cstdint>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/create_basic.hpp"

using namespace boost::mysql::detail;
using namespace boost::unit_test;
using namespace boost::mysql::test;
using boost::mysql::field_view;

BOOST_AUTO_TEST_SUITE(test_null_bitmap)

BOOST_AUTO_TEST_SUITE(parser)

BOOST_AUTO_TEST_CASE(byte_count)
{
    struct
    {
        std::size_t num_fields;
        std::size_t expected;
    } test_cases[] = {
        {0,  1},
        {1,  1},
        {2,  1},
        {3,  1},
        {4,  1},
        {5,  1},
        {6,  1},
        {7,  2},
        {8,  2},
        {9,  2},
        {10, 2},
        {11, 2},
        {12, 2},
        {13, 2},
        {14, 2},
        {15, 3},
        {16, 3},
        {17, 3},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.num_fields)
        {
            null_bitmap_parser parser{tc.num_fields};
            BOOST_TEST(parser.byte_count() == tc.expected);
        }
    }
}

BOOST_AUTO_TEST_CASE(is_null_coverage)
{
    // Given a null bitmap with 17 fields, and the following buffer, vary the field offset
    // 0b10110100, 0b11111111, 0b00000000
    struct
    {
        std::size_t pos;
        bool expected;
    } test_cases[] = {
        {0,  true },
        {1,  false},
        {2,  true },
        {3,  true },
        {4,  false},
        {5,  true },
        {6,  true },
        {7,  true },
        {8,  true },
        {9,  true },
        {10, true },
        {11, true },
        {12, true },
        {13, true },
        {14, false},
        {15, false},
        {16, false},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.pos)
        {
            const std::array<std::uint8_t, 3> buffer{
                {0xb4, 0xff, 0x00}
            };
            null_bitmap_parser parser{17};  // 17 fields
            bool actual = parser.is_null(buffer.data(), tc.pos);
            BOOST_TEST(actual == tc.expected);
        }
    }
}

// spotcheck: we handle the offset correctly, ignoring the first two bits
BOOST_AUTO_TEST_CASE(is_null_first_bit)
{
    struct
    {
        std::uint8_t buffer;
        bool expected;
    } test_cases[] = {
        {0x00, false},
        {0x01, false},
        {0x02, false},
        {0x03, false},
        {0x04, true },
        {0x05, true },
        {0x06, true },
        {0x07, true },
        {0x08, false},
        {0x09, false},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.buffer)
        {
            null_bitmap_parser parser{1};  // 1 field
            bool actual = parser.is_null(&tc.buffer, 0);
            BOOST_TEST(actual == tc.expected);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(generator)

BOOST_AUTO_TEST_CASE(coverage)
{
    // Creates a generator and calls it until exhausted
    auto gen_null_bitmap = [](boost::span<const field_view> params) {
        std::vector<std::uint8_t> res;
        null_bitmap_generator gen(params);
        while (!gen.done())
            res.push_back(gen.next());
        return res;
    };

    struct
    {
        const char* name;
        std::vector<field_view> params;
        std::vector<std::uint8_t> expected;
    } test_cases[] = {
        // All combinations for up to 2 values
        {"empty", {}, {}},
        {"N", make_fv_vector(nullptr), {0x01}},
        {"V", make_fv_vector(1), {0x00}},
        {"NV", make_fv_vector(nullptr, 1), {0x01}},
        {"NN", make_fv_vector(nullptr, nullptr), {0x03}},
        {"VV", make_fv_vector(1, 1), {0x00}},
        {"VN", make_fv_vector(1, nullptr), {0x02}},

        // Last value null - checking we set the right bit
        {"VVN", make_fv_vector(1, 1, nullptr), {0x04}},
        {"VVVN", make_fv_vector(1, 1, 1, nullptr), {0x08}},
        {"VVVVN", make_fv_vector(1, 1, 1, 1, nullptr), {0x10}},
        {"VVVVVN", make_fv_vector(1, 1, 1, 1, 1, nullptr), {0x20}},
        {"VVVVVVN", make_fv_vector(1, 1, 1, 1, 1, 1, nullptr), {0x40}},
        {"VVVVVVVN", make_fv_vector(1, 1, 1, 1, 1, 1, 1, nullptr), {0x80}},
        {"VVVVVVVV_N", make_fv_vector(1, 1, 1, 1, 1, 1, 1, 1, nullptr), {0x00, 0x01}},
        {"VVVVVVVV_VN", make_fv_vector(1, 1, 1, 1, 1, 1, 1, 1, 1, nullptr), {0x00, 0x02}},
        {"VVVVVVVV_VVN", make_fv_vector(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, nullptr), {0x00, 0x04}},
        {"VVVVVVVV_VVVN", make_fv_vector(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, nullptr), {0x00, 0x08}},
        {"VVVVVVVV_VVVVN", make_fv_vector(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, nullptr), {0x00, 0x10}},
        {"VVVVVVVV_VVVVVN", make_fv_vector(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, nullptr), {0x00, 0x20}},
        {"VVVVVVVV_VVVVVVN", make_fv_vector(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, nullptr), {0x00, 0x40}},
        {"VVVVVVVV_VVVVVVVN",
         make_fv_vector(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, nullptr),
         {0x00, 0x80}},
        {"VVVVVVVV_VVVVVVVV_N",
         make_fv_vector(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, nullptr),
         {0x00, 0x00, 0x01}},

        // Some arbitrary combinations
        {"VNVVNVVN", make_fv_vector(1, nullptr, 1, 1, nullptr, 1, 1, nullptr), {0x92}},
        {"NVVVVNVV_VVNVN",
         make_fv_vector(nullptr, 1, 1, 1, 1, nullptr, 1, 1, 1, 1, nullptr, 1, nullptr),
         {0x21, 0x14}},
        {"VVVVVVVV_VVVVVVVV_V",
         make_fv_vector(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
         {0x00, 0x00, 0x00}},
        {"NNNNNNNN_NNNNNNNN_NNN", std::vector<field_view>(19, field_view()), {0xff, 0xff, 0x07}},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            auto actual = gen_null_bitmap(tc.params);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(actual, tc.expected);
        }
    }
}

// Spotcheck: generating step-by-step yields the expected results
BOOST_AUTO_TEST_CASE(step_by_step)
{
    // Setup
    auto params = make_fv_arr(
        // byte 1
        nullptr,
        1,
        1,
        1,
        nullptr,
        nullptr,
        1,
        1,

        // byte 2
        1,
        nullptr,
        nullptr,
        1,
        1,
        1,
        1,
        nullptr,

        // byte 3
        1,
        1,
        1,
        nullptr,
        1,
        nullptr
    );
    null_bitmap_generator gen(params);

    // Initiates as not done
    BOOST_TEST(!gen.done());

    // Generate first byte
    BOOST_TEST(gen.next() == 0x31u);
    BOOST_TEST(!gen.done());

    // Generate second byte
    BOOST_TEST(gen.next() == 0x86u);
    BOOST_TEST(!gen.done());

    // Generate last byte
    BOOST_TEST(gen.next() == 0x28u);
    BOOST_TEST(gen.done());
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()  // test_null_bitmap_traits
