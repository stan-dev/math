//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_DETAIL_PROTOCOL_SERIALIZATION_TEST_SAMPLES_BASIC_TYPES_HPP
#define BOOST_MYSQL_TEST_UNIT_DETAIL_PROTOCOL_SERIALIZATION_TEST_SAMPLES_BASIC_TYPES_HPP

#include <boost/mysql/detail/protocol/serialization.hpp>

#include "../serialization_test.hpp"
#include "buffer_concat.hpp"

namespace boost {
namespace mysql {
namespace test {

// Definitions for the parameterized tests
const std::string string_250(250, 'a');
const std::string string_251(251, 'a');
const std::string string_ffff(0xffff, 'a');
const std::string string_10000(0x10000, 'a');

enum class enum_int1 : std::uint8_t
{
    value0 = 0,
    value1 = 3,
    value2 = 0xff
};

enum class enum_int2 : std::uint16_t
{
    value0 = 0,
    value1 = 3,
    value2 = 0xfeff
};

enum class enum_int4 : std::uint32_t
{
    value0 = 0,
    value1 = 3,
    value2 = 0xfcfdfeff
};

// clang-format off
const serialization_test_spec int_spec {
    serialization_test_type::full, {
        { "int1", std::uint8_t(0xff), {0xff} },
        { "int2", std::uint16_t(0xfeff), {0xff, 0xfe} },
        { "int3", int3(0xfdfeff), {0xff, 0xfe, 0xfd} },
        { "int4", std::uint32_t(0xfcfdfeff), {0xff, 0xfe, 0xfd, 0xfc} },
        { "int8", std::uint64_t(0xf8f9fafbfcfdfeff), {0xff, 0xfe, 0xfd, 0xfc, 0xfb, 0xfa, 0xf9, 0xf8} },
        { "int1_positive", std::int8_t(0x01), {0x01} },
        { "int1_negative", std::int8_t(-1), {0xff} },
        { "int2_positive", std::int16_t(0x0201), {0x01, 0x02} },
        { "int2_negative", std::int16_t(-0x101), {0xff, 0xfe} },
        { "int4_positive", std::int32_t(0x04030201), {0x01, 0x02, 0x03, 0x04} },
        { "int4_negative", std::int32_t(-0x3020101), {0xff, 0xfe, 0xfd, 0xfc} },
        { "int8_positive", std::int64_t(0x0807060504030201),
                {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08} },
        { "int8_negative", std::int64_t(-0x0706050403020101),
                {0xff, 0xfe, 0xfd, 0xfc, 0xfb, 0xfa, 0xf9, 0xf8} },

        // int_lenenc
        { "1_byte_regular", int_lenenc(1), {0x01} },
        { "1_byte_max", int_lenenc(250), {0xfa} },
        { "2_bytes_regular", int_lenenc(0xfeb7), {0xfc, 0xb7, 0xfe} },
        { "2_bytes_max", int_lenenc(0xffff), {0xfc, 0xff, 0xff} },
        { "3_bytes_regular", int_lenenc(0xa0feff), {0xfd, 0xff, 0xfe, 0xa0} },
        { "3_bytes_max", int_lenenc(0xffffff), {0xfd, 0xff, 0xff, 0xff} },
        { "8_bytes_regular", int_lenenc(0xf8f9fafbfcfdfeff),
                {0xfe, 0xff, 0xfe, 0xfd, 0xfc, 0xfb, 0xfa, 0xf9, 0xf8} },
        { "8_bytes_max", int_lenenc(0xffffffffffffffff),
                {0xfe, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff} }
    }
};

const serialization_test_spec enum_spec {
    serialization_test_type::full, {
        { "enum_int1_low_value", enum_int1::value1, {0x03} },
        { "enum_int1_high_value", enum_int1::value2, {0xff} },
        { "enum_int2_low_value", enum_int2::value1, {0x03, 0x00} },
        { "enum_int2_high_value", enum_int2::value2, {0xff, 0xfe} },
        { "enum_int4_low_value", enum_int4::value1, {0x03, 0x00, 0x00, 0x00} },
        { "enum_int4_high_value", enum_int4::value2, {0xff, 0xfe, 0xfd, 0xfc} }
    }
};

const serialization_test_spec string_fixed_spec {
    serialization_test_type::full, {
        { "4c_regular_characters", makesfixed<4>("abde"), {0x61, 0x62, 0x64, 0x65} },
        { "3c_null_characters", makesfixed<3>("\0\1a"), {0x00, 0x01, 0x61} },
        { "3c_utf8_characters", string_fixed<3>{{char(0xc3), char(0xb1), 'a'}}, {0xc3, 0xb1, 0x61} },
        { "1c_regular_characters", makesfixed<1>("a"), {0x61} }
    }
};

const serialization_test_spec string_null_spec {
    serialization_test_type::full, {
        { "regular_characters", string_null("abc"), {0x61, 0x62, 0x63, 0x00} },
        { "utf8_characters", string_null("\xc3\xb1"), {0xc3, 0xb1, 0x00} },
        { "empty", string_null(""), {0x00} }
    }
};

const serialization_test_spec string_lenenc_spec {
    serialization_test_type::full, {
        { "empty", string_lenenc(""),
            {0x00} },
        { "1_byte_size_regular_characters", string_lenenc("abc"),
            {0x03, 0x61, 0x62, 0x63} },
        { "1_byte_size_null_characters", string_lenenc(makesv("a\0b")),
            {0x03, 0x61, 0x00, 0x62} },
        { "1_byte_size_max", string_lenenc(string_250),
            concat_copy({250}, std::vector<std::uint8_t>(250, 0x61)) },
        { "2_byte_size_min", string_lenenc(string_251),
            concat_copy({0xfc, 251, 0}, std::vector<std::uint8_t>(251, 0x61)) },
        { "2_byte_size_max", string_lenenc(string_ffff),
            concat_copy({0xfc, 0xff, 0xff}, std::vector<std::uint8_t>(0xffff, 0x61)) },
        { "3_byte_size_min", string_lenenc(string_10000),
            concat_copy({0xfd, 0x00, 0x00, 0x01}, std::vector<std::uint8_t>(0x10000, 0x61)) }
    }
};

const serialization_test_spec string_eof_spec {
    serialization_test_type::full_no_space, {
        { "regular_characters", string_eof("abc"), {0x61, 0x62, 0x63} },
        { "null_characters", string_eof(string_view("a\0b", 3)), {0x61, 0x00, 0x62} },
        { "empty", string_eof(""), {} }
    }
};
// clang-format on

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
