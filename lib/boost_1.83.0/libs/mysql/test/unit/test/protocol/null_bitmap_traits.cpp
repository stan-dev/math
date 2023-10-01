//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/impl/internal/protocol/null_bitmap_traits.hpp>

#include <boost/test/data/monomorphic/collection.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <array>

using namespace boost::mysql::detail;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(test_null_bitmap_traits)

// byte_count
struct byte_count_sample
{
    std::size_t offset;
    std::size_t num_fields;
    std::size_t expected_value;
};

std::ostream& operator<<(std::ostream& os, const byte_count_sample& val)
{
    return os << "(offset=" << val.offset << ", num_fields=" << val.num_fields << ")";
}

constexpr byte_count_sample all_byte_count_samples[]{
    {stmt_execute_null_bitmap_offset, 0,  0},
    {stmt_execute_null_bitmap_offset, 1,  1},
    {stmt_execute_null_bitmap_offset, 2,  1},
    {stmt_execute_null_bitmap_offset, 3,  1},
    {stmt_execute_null_bitmap_offset, 4,  1},
    {stmt_execute_null_bitmap_offset, 5,  1},
    {stmt_execute_null_bitmap_offset, 6,  1},
    {stmt_execute_null_bitmap_offset, 7,  1},
    {stmt_execute_null_bitmap_offset, 8,  1},
    {stmt_execute_null_bitmap_offset, 9,  2},
    {stmt_execute_null_bitmap_offset, 10, 2},
    {stmt_execute_null_bitmap_offset, 11, 2},
    {stmt_execute_null_bitmap_offset, 12, 2},
    {stmt_execute_null_bitmap_offset, 13, 2},
    {stmt_execute_null_bitmap_offset, 14, 2},
    {stmt_execute_null_bitmap_offset, 15, 2},
    {stmt_execute_null_bitmap_offset, 16, 2},
    {stmt_execute_null_bitmap_offset, 17, 3},

    {binary_row_null_bitmap_offset,   0,  1},
    {binary_row_null_bitmap_offset,   1,  1},
    {binary_row_null_bitmap_offset,   2,  1},
    {binary_row_null_bitmap_offset,   3,  1},
    {binary_row_null_bitmap_offset,   4,  1},
    {binary_row_null_bitmap_offset,   5,  1},
    {binary_row_null_bitmap_offset,   6,  1},
    {binary_row_null_bitmap_offset,   7,  2},
    {binary_row_null_bitmap_offset,   8,  2},
    {binary_row_null_bitmap_offset,   9,  2},
    {binary_row_null_bitmap_offset,   10, 2},
    {binary_row_null_bitmap_offset,   11, 2},
    {binary_row_null_bitmap_offset,   12, 2},
    {binary_row_null_bitmap_offset,   13, 2},
    {binary_row_null_bitmap_offset,   14, 2},
    {binary_row_null_bitmap_offset,   15, 3},
    {binary_row_null_bitmap_offset,   16, 3},
    {binary_row_null_bitmap_offset,   17, 3},
};

BOOST_DATA_TEST_CASE(byte_count, data::make(all_byte_count_samples))
{
    null_bitmap_traits traits(sample.offset, sample.num_fields);
    BOOST_TEST(traits.byte_count() == sample.expected_value);
}

// is_null
struct is_null_sample
{
    std::size_t offset;
    std::size_t pos;
    bool expected;
};

std::ostream& operator<<(std::ostream& os, const is_null_sample& value)
{
    return os << "(offset=" << value.offset << ", pos=" << value.pos << ")";
}

constexpr is_null_sample all_is_null_samples[]{
    {stmt_execute_null_bitmap_offset, 0,  false},
    {stmt_execute_null_bitmap_offset, 1,  false},
    {stmt_execute_null_bitmap_offset, 2,  true },
    {stmt_execute_null_bitmap_offset, 3,  false},
    {stmt_execute_null_bitmap_offset, 4,  true },
    {stmt_execute_null_bitmap_offset, 5,  true },
    {stmt_execute_null_bitmap_offset, 6,  false},
    {stmt_execute_null_bitmap_offset, 7,  true },
    {stmt_execute_null_bitmap_offset, 8,  true },
    {stmt_execute_null_bitmap_offset, 9,  true },
    {stmt_execute_null_bitmap_offset, 10, true },
    {stmt_execute_null_bitmap_offset, 11, true },
    {stmt_execute_null_bitmap_offset, 12, true },
    {stmt_execute_null_bitmap_offset, 13, true },
    {stmt_execute_null_bitmap_offset, 14, true },
    {stmt_execute_null_bitmap_offset, 15, true },
    {stmt_execute_null_bitmap_offset, 16, false},

    {binary_row_null_bitmap_offset,   0,  true },
    {binary_row_null_bitmap_offset,   1,  false},
    {binary_row_null_bitmap_offset,   2,  true },
    {binary_row_null_bitmap_offset,   3,  true },
    {binary_row_null_bitmap_offset,   4,  false},
    {binary_row_null_bitmap_offset,   5,  true },
    {binary_row_null_bitmap_offset,   6,  true },
    {binary_row_null_bitmap_offset,   7,  true },
    {binary_row_null_bitmap_offset,   8,  true },
    {binary_row_null_bitmap_offset,   9,  true },
    {binary_row_null_bitmap_offset,   10, true },
    {binary_row_null_bitmap_offset,   11, true },
    {binary_row_null_bitmap_offset,   12, true },
    {binary_row_null_bitmap_offset,   13, true },
    {binary_row_null_bitmap_offset,   14, false},
    {binary_row_null_bitmap_offset,   15, false},
    {binary_row_null_bitmap_offset,   16, false},
};

BOOST_DATA_TEST_CASE(is_null, data::make(all_is_null_samples))
{
    // 0b10110100, 0b11111111, 0b00000000
    std::uint8_t content[] = {0xb4, 0xff, 0x00};
    null_bitmap_traits traits(sample.offset, 17);  // 17 fields
    bool actual = traits.is_null(content, sample.pos);
    BOOST_TEST(actual == sample.expected);
}

BOOST_AUTO_TEST_CASE(is_null_one_field_stmt_execute_first_bit_zero)
{
    std::uint8_t value = 0x00;
    null_bitmap_traits traits(stmt_execute_null_bitmap_offset, 1);
    BOOST_TEST(!traits.is_null(&value, 0));
}

BOOST_AUTO_TEST_CASE(is_null_one_field_stmt_execute_first_bit_one)
{
    std::uint8_t value = 0x01;
    null_bitmap_traits traits(stmt_execute_null_bitmap_offset, 1);
    BOOST_TEST(traits.is_null(&value, 0));
}

BOOST_AUTO_TEST_CASE(is_null_one_field_binary_row_third_bit_zero)
{
    std::uint8_t value = 0x00;
    null_bitmap_traits traits(binary_row_null_bitmap_offset, 1);
    BOOST_TEST(!traits.is_null(&value, 0));
}

BOOST_AUTO_TEST_CASE(is_null_one_field_binary_row_third_bit_one)
{
    std::uint8_t value = 0x04;  // 0b00000100
    null_bitmap_traits traits(binary_row_null_bitmap_offset, 1);
    BOOST_TEST(traits.is_null(&value, 0));
}

// set_null
struct set_null_sample
{
    std::size_t offset;
    std::size_t pos;
    std::array<std::uint8_t, 3> expected;

    constexpr set_null_sample(
        std::size_t offset,
        std::size_t pos,
        const std::array<std::uint8_t, 3>& expected
    )
        : offset(offset), pos(pos), expected(expected){};
};

std::ostream& operator<<(std::ostream& os, const set_null_sample& value)
{
    return os << "(offset=" << value.offset << ", pos=" << value.pos << ")";
}

constexpr set_null_sample all_set_null_samples[]{
    {stmt_execute_null_bitmap_offset, 0,  {{0x01, 0, 0}}},
    {stmt_execute_null_bitmap_offset, 1,  {{0x02, 0, 0}}},
    {stmt_execute_null_bitmap_offset, 2,  {{0x04, 0, 0}}},
    {stmt_execute_null_bitmap_offset, 3,  {{0x08, 0, 0}}},
    {stmt_execute_null_bitmap_offset, 4,  {{0x10, 0, 0}}},
    {stmt_execute_null_bitmap_offset, 5,  {{0x20, 0, 0}}},
    {stmt_execute_null_bitmap_offset, 6,  {{0x40, 0, 0}}},
    {stmt_execute_null_bitmap_offset, 7,  {{0x80, 0, 0}}},
    {stmt_execute_null_bitmap_offset, 8,  {{0, 0x01, 0}}},
    {stmt_execute_null_bitmap_offset, 9,  {{0, 0x02, 0}}},
    {stmt_execute_null_bitmap_offset, 10, {{0, 0x04, 0}}},
    {stmt_execute_null_bitmap_offset, 11, {{0, 0x08, 0}}},
    {stmt_execute_null_bitmap_offset, 12, {{0, 0x10, 0}}},
    {stmt_execute_null_bitmap_offset, 13, {{0, 0x20, 0}}},
    {stmt_execute_null_bitmap_offset, 14, {{0, 0x40, 0}}},
    {stmt_execute_null_bitmap_offset, 15, {{0, 0x80, 0}}},
    {stmt_execute_null_bitmap_offset, 16, {{0, 0, 0x01}}},

    {binary_row_null_bitmap_offset,   0,  {{0x04, 0, 0}}},
    {binary_row_null_bitmap_offset,   1,  {{0x08, 0, 0}}},
    {binary_row_null_bitmap_offset,   2,  {{0x10, 0, 0}}},
    {binary_row_null_bitmap_offset,   3,  {{0x20, 0, 0}}},
    {binary_row_null_bitmap_offset,   4,  {{0x40, 0, 0}}},
    {binary_row_null_bitmap_offset,   5,  {{0x80, 0, 0}}},
    {binary_row_null_bitmap_offset,   6,  {{0, 0x01, 0}}},
    {binary_row_null_bitmap_offset,   7,  {{0, 0x02, 0}}},
    {binary_row_null_bitmap_offset,   8,  {{0, 0x04, 0}}},
    {binary_row_null_bitmap_offset,   9,  {{0, 0x08, 0}}},
    {binary_row_null_bitmap_offset,   10, {{0, 0x10, 0}}},
    {binary_row_null_bitmap_offset,   11, {{0, 0x20, 0}}},
    {binary_row_null_bitmap_offset,   12, {{0, 0x40, 0}}},
    {binary_row_null_bitmap_offset,   13, {{0, 0x80, 0}}},
    {binary_row_null_bitmap_offset,   14, {{0, 0, 0x01}}},
    {binary_row_null_bitmap_offset,   15, {{0, 0, 0x02}}},
    {binary_row_null_bitmap_offset,   16, {{0, 0, 0x04}}},
};

BOOST_DATA_TEST_CASE(set_null, data::make(all_set_null_samples))
{
    std::array<std::uint8_t, 4> expected_buffer{};  // help detect buffer overruns
    std::memcpy(expected_buffer.data(), sample.expected.data(), 3);
    std::array<std::uint8_t, 4> actual_buffer{};
    null_bitmap_traits traits(sample.offset, 17);  // 17 fields
    traits.set_null(actual_buffer.data(), sample.pos);
    BOOST_TEST(expected_buffer == actual_buffer);
}

BOOST_AUTO_TEST_CASE(set_null_one_field_stmt_execute)
{
    std::uint8_t value = 0;
    null_bitmap_traits traits(stmt_execute_null_bitmap_offset, 1);
    traits.set_null(&value, 0);
    BOOST_TEST(value == 1);
}

BOOST_AUTO_TEST_CASE(set_null_one_field_binary_row)
{
    std::uint8_t value = 0;
    null_bitmap_traits traits(binary_row_null_bitmap_offset, 1);
    traits.set_null(&value, 0);
    BOOST_TEST(value == 4);
}

BOOST_AUTO_TEST_CASE(set_null_multified_stmt_execute)
{
    std::array<std::uint8_t, 4> expected_buffer{
        {0xb4, 0xff, 0x00, 0x00}
    };
    std::array<std::uint8_t, 4> actual_buffer{};
    null_bitmap_traits traits(stmt_execute_null_bitmap_offset, 17);  // 17 fields
    for (std::size_t pos : {2, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15})
    {
        traits.set_null(actual_buffer.data(), pos);
    }
    BOOST_TEST(expected_buffer == actual_buffer);
}

BOOST_AUTO_TEST_CASE(set_null_multified_binary_row)
{
    std::array<std::uint8_t, 4> expected_buffer{
        {0xb4, 0xff, 0x00, 0x00}
    };
    std::array<std::uint8_t, 4> actual_buffer{};
    null_bitmap_traits traits(binary_row_null_bitmap_offset, 17);  // 17 fields
    for (std::size_t pos : {0, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13})
    {
        traits.set_null(actual_buffer.data(), pos);
    }
    BOOST_TEST(expected_buffer == actual_buffer);
}

BOOST_AUTO_TEST_SUITE_END()  // test_null_bitmap_traits
