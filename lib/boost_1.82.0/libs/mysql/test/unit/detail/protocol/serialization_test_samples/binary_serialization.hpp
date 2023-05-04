//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_DETAIL_PROTOCOL_SERIALIZATION_TEST_SAMPLES_BINARY_SERIALIZATION_HPP
#define BOOST_MYSQL_TEST_UNIT_DETAIL_PROTOCOL_SERIALIZATION_TEST_SAMPLES_BINARY_SERIALIZATION_HPP

#include <boost/mysql/blob_view.hpp>

#include <boost/mysql/detail/protocol/binary_serialization.hpp>

#include "../serialization_test.hpp"

namespace boost {
namespace mysql {
namespace test {

template <class T>
serialization_sample make_binary_serialization_sample(
    std::string&& name,
    const T& expected,
    detail::bytestring&& buffer
)
{
    return serialization_sample(std::move(name), field_view(expected), std::move(buffer));
}

constexpr std::uint8_t binary_serialization_blob_buffer[] = {0x70, 0x00, 0x01};

// clang-format off
const serialization_test_spec binary_serialization_spec {
    serialization_test_type::serialization, {
        // Strings and ints: extensive testing already done. Ensure
        // we call the right function
        make_binary_serialization_sample("string", "abc",
            {0x03, 0x61, 0x62, 0x63}),
        make_binary_serialization_sample("blob", blob_view(binary_serialization_blob_buffer),
            {0x03, 0x70, 0x00, 0x01}),
        make_binary_serialization_sample("uint64", std::uint64_t(0xf8f9fafbfcfdfeff),
            {0xff, 0xfe, 0xfd, 0xfc, 0xfb, 0xfa, 0xf9, 0xf8}),
        make_binary_serialization_sample("int64",  std::int64_t(-0x0706050403020101),
            {0xff, 0xfe, 0xfd, 0xfc, 0xfb, 0xfa, 0xf9, 0xf8}),
        make_binary_serialization_sample("float_fractional_negative", -4.2f,
            {0x66, 0x66, 0x86, 0xc0}),
        make_binary_serialization_sample("float_fractional_positive", 4.2f,
            {0x66, 0x66, 0x86, 0x40}),
        make_binary_serialization_sample("float_positive_exp_positive_fractional", 3.14e20f,
            {0x01, 0x2d, 0x88, 0x61}),
        make_binary_serialization_sample("float_zero", 0.0f,
            {0x00, 0x00, 0x00, 0x00}),
        make_binary_serialization_sample("double_fractional_negative", -4.2,
            {0xcd, 0xcc, 0xcc, 0xcc, 0xcc, 0xcc, 0x10, 0xc0}),
        make_binary_serialization_sample("double_fractional_positive", 4.2,
            {0xcd, 0xcc, 0xcc, 0xcc, 0xcc, 0xcc, 0x10, 0x40}),
        make_binary_serialization_sample("double_positive_exp_positive_fractional", 3.14e200,
            {0xce, 0x46, 0x3c, 0x76, 0x9c, 0x68, 0x90, 0x69}),
        make_binary_serialization_sample("double_zero", 0.0,
            {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}),
        make_binary_serialization_sample("date_regular", date(2010u, 3u, 28u),
            {0x04, 0xda, 0x07, 0x03, 0x1c}),
        make_binary_serialization_sample("date_min", date(1000u, 1u, 1u),
            {0x04, 0xe8, 0x03, 0x01, 0x01}),
        make_binary_serialization_sample("date_mysqlmax", date(9999u, 12u, 31u),
            {0x04, 0x0f, 0x27, 0x0c, 0x1f}),
        make_binary_serialization_sample("date_max", date(0xffff, 0xff, 0xff),
            {0x04, 0xff, 0xff, 0xff, 0xff}),
        make_binary_serialization_sample("date_zero", date(),
            {0x04, 0x00, 0x00, 0x00, 0x00}),
        make_binary_serialization_sample("datetime_regular", datetime(2010u, 1u, 1u, 23u,  1u,  59u, 967510u),
            {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00}),
        make_binary_serialization_sample("datetime_min", datetime(0, 1, 1),
            {0x0b, 0x00, 0x00, 0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}),
        make_binary_serialization_sample("datetime_mysqlmax", datetime(9999, 12, 31, 23, 59, 59, 999999),
            {0x0b, 0x0f, 0x27, 0x0c, 0x1f, 0x17, 0x3b, 0x3b, 0x3f, 0x42, 0x0f, 0x00}),
        make_binary_serialization_sample("datetime_max", datetime(0xffff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xffffffff),
            {0x0b, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff}),
        make_binary_serialization_sample("datetime_zero", datetime(),
            {0x0b, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}),
        make_binary_serialization_sample("time_positive_u", maket(0, 0, 0, 321000),
            {0x0c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xe8, 0xe5, 0x04, 0x00}),
        make_binary_serialization_sample("time_positive_hmsu", maket(838, 59, 58, 999000),
            {0x0c, 0x00, 0x22, 0x00, 0x00, 0x00, 0x16, 0x3b, 0x3a, 0x58, 0x3e, 0x0f, 0x00}),
        make_binary_serialization_sample("time_negative_u", -maket(0, 0, 0, 321000),
            {0x0c, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xe8, 0xe5, 0x04, 0x00}),
        make_binary_serialization_sample("time_negative_hmsu", -maket(838, 59, 58, 999000),
            {0x0c, 0x01, 0x22, 0x00, 0x00, 0x00, 0x16, 0x3b, 0x3a, 0x58, 0x3e, 0x0f, 0x00}),

        // NULL is transmitted as the NULL bitmap, so nothing is expected as output
        make_binary_serialization_sample("null", nullptr, {}),
    }
};
// clang-format on

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
