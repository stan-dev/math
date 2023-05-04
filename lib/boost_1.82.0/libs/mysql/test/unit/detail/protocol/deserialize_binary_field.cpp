//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/blob_view.hpp>
#include <boost/mysql/date.hpp>
#include <boost/mysql/datetime.hpp>
#include <boost/mysql/mysql_collations.hpp>

#include <boost/mysql/detail/auxiliar/access_fwd.hpp>
#include <boost/mysql/detail/auxiliar/stringize.hpp>
#include <boost/mysql/detail/protocol/constants.hpp>
#include <boost/mysql/detail/protocol/deserialize_binary_field.hpp>
#include <boost/mysql/detail/protocol/deserialize_errc.hpp>

#include <boost/test/data/monomorphic/collection.hpp>
#include <boost/test/data/test_case.hpp>

#include <cstddef>

#include "create_meta.hpp"
#include "printing.hpp"
#include "test_common.hpp"

using namespace boost::mysql::detail;
using namespace boost::mysql::test;
using namespace boost::unit_test;
using boost::mysql::blob_view;
using boost::mysql::date;
using boost::mysql::datetime;
using boost::mysql::error_code;
using boost::mysql::field_view;

namespace {

BOOST_AUTO_TEST_SUITE(test_deserialize_binary_field)

BOOST_AUTO_TEST_SUITE(success)

struct success_sample
{
    std::string name;
    std::vector<std::uint8_t> from;
    field_view expected;
    protocol_field_type type;
    std::uint16_t flags;
    std::uint16_t collation;

    template <class T>
    success_sample(
        std::string name,
        std::vector<std::uint8_t> from,
        T&& expected_value,
        protocol_field_type type,
        std::uint16_t flags = 0,
        std::uint16_t collation = boost::mysql::mysql_collations::utf8mb4_general_ci
    )
        : name(std::move(name)),
          from(std::move(from)),
          expected(std::forward<T>(expected_value)),
          type(type),
          flags(flags),
          collation(collation)
    {
    }
};

std::ostream& operator<<(std::ostream& os, const success_sample& input)
{
    return os << "(type=" << input.type << ", name=" << input.name << ")";
}

void add_string_samples(std::vector<success_sample>& output)
{
    output.push_back(
        success_sample("varchar", {0x04, 0x74, 0x65, 0x73, 0x74}, "test", protocol_field_type::var_string)
    );
    output.push_back(
        success_sample("char", {0x04, 0x74, 0x65, 0x73, 0x74}, "test", protocol_field_type::string)
    );
    output.push_back(success_sample(
        "text",
        {0x04, 0x74, 0x65, 0x73, 0x74},
        "test",
        protocol_field_type::blob,
        column_flags::blob
    ));
    output.push_back(success_sample(
        "enum",
        {0x04, 0x74, 0x65, 0x73, 0x74},
        "test",
        protocol_field_type::string,
        column_flags::enum_
    ));
    output.push_back(success_sample(
        "set",
        {0x04, 0x74, 0x65, 0x73, 0x74},
        "test",
        protocol_field_type::string,
        column_flags::set
    ));
    output.push_back(success_sample("decimal", {0x02, 0x31, 0x30}, "10", protocol_field_type::newdecimal));
    output.push_back(success_sample("json", {0x02, 0x7b, 0x7d}, "{}", protocol_field_type::json));
}

void add_blob_samples(std::vector<success_sample>& output)
{
    static constexpr std::uint8_t buff[] = {0x01, 0x00, 0x73, 0x74};

    output.push_back(success_sample(
        "varbinary",
        {0x04, 0x01, 0x00, 0x73, 0x74},
        blob_view(buff),
        protocol_field_type::var_string,
        column_flags::binary,
        binary_collation
    ));
    output.push_back(success_sample(
        "binary",
        {0x04, 0x01, 0x00, 0x73, 0x74},
        blob_view(buff),
        protocol_field_type::string,
        column_flags::binary,
        binary_collation
    ));
    output.push_back(success_sample(
        "blob",
        {0x04, 0x01, 0x00, 0x73, 0x74},
        blob_view(buff),
        protocol_field_type::blob,
        column_flags::binary,
        binary_collation
    ));
    output.push_back(success_sample(
        "geometry",
        {0x04, 0x01, 0x00, 0x73, 0x74},
        blob_view(buff),
        protocol_field_type::geometry,
        binary_collation
    ));

    // Anything we don't know what it is, we interpret as a blob
    output.push_back(success_sample(
        "unknown_protocol_type",
        {0x04, 0x01, 0x00, 0x73, 0x74},
        blob_view(buff),
        static_cast<protocol_field_type>(0x23)
    ));
}

// Note: these employ regular integer deserialization functions, which have
// already been tested
void add_int_samples(std::vector<success_sample>& output)
{
    output.push_back(success_sample(
        "tinyint_unsigned",
        {0x14},
        std::uint64_t(20),
        protocol_field_type::tiny,
        column_flags::unsigned_
    ));
    output.push_back(success_sample("tinyint_signed", {0xec}, std::int64_t(-20), protocol_field_type::tiny));

    output.push_back(success_sample(
        "smallint_unsigned",
        {0x14, 0x00},
        std::uint64_t(20),
        protocol_field_type::short_,
        column_flags::unsigned_
    ));
    output.push_back(
        success_sample("smallint_signed", {0xec, 0xff}, std::int64_t(-20), protocol_field_type::short_)
    );

    output.push_back(success_sample(
        "mediumint_unsigned",
        {0x14, 0x00, 0x00, 0x00},
        std::uint64_t(20),
        protocol_field_type::int24,
        column_flags::unsigned_
    ));
    output.push_back(success_sample(
        "mediumint_signed",
        {0xec, 0xff, 0xff, 0xff},
        std::int64_t(-20),
        protocol_field_type::int24
    ));

    output.push_back(success_sample(
        "int_unsigned",
        {0x14, 0x00, 0x00, 0x00},
        std::uint64_t(20),
        protocol_field_type::long_,
        column_flags::unsigned_
    ));
    output.push_back(
        success_sample("int_signed", {0xec, 0xff, 0xff, 0xff}, std::int64_t(-20), protocol_field_type::long_)
    );

    output.push_back(success_sample(
        "bigint_unsigned",
        {0x14, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
        std::uint64_t(20),
        protocol_field_type::longlong,
        column_flags::unsigned_
    ));
    output.push_back(success_sample(
        "bigint_signed",
        {0xec, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff},
        std::int64_t(-20),
        protocol_field_type::longlong
    ));

    output.push_back(success_sample(
        "year",
        {0xe3, 0x07},
        std::uint64_t(2019),
        protocol_field_type::year,
        column_flags::unsigned_
    ));
}

// bit
void add_bit_types(std::vector<success_sample>& output)
{
    output.push_back(success_sample(
        "bit_8",
        {0x01, 0x12},
        std::uint64_t(0x12),
        protocol_field_type::bit,
        column_flags::unsigned_
    ));
    output.push_back(success_sample(
        "bit_16",
        {0x02, 0x12, 0x34},
        std::uint64_t(0x1234),
        protocol_field_type::bit,
        column_flags::unsigned_
    ));
    output.push_back(success_sample(
        "bit_24",
        {0x03, 0x12, 0x34, 0x56},
        std::uint64_t(0x123456),
        protocol_field_type::bit,
        column_flags::unsigned_
    ));
    output.push_back(success_sample(
        "bit_32",
        {0x04, 0x12, 0x34, 0x56, 0x78},
        std::uint64_t(0x12345678),
        protocol_field_type::bit,
        column_flags::unsigned_
    ));
    output.push_back(success_sample(
        "bit_40",
        {0x05, 0x12, 0x34, 0x56, 0x78, 0x9a},
        std::uint64_t(0x123456789a),
        protocol_field_type::bit,
        column_flags::unsigned_
    ));
    output.push_back(success_sample(
        "bit_48",
        {0x06, 0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc},
        std::uint64_t(0x123456789abc),
        protocol_field_type::bit,
        column_flags::unsigned_
    ));
    output.push_back(success_sample(
        "bit_56",
        {0x07, 0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde},
        std::uint64_t(0x123456789abcde),
        protocol_field_type::bit,
        column_flags::unsigned_
    ));
    output.push_back(success_sample(
        "bit_64",
        {0x08, 0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xf0},
        std::uint64_t(0x123456789abcdef0),
        protocol_field_type::bit,
        column_flags::unsigned_
    ));
}

void add_float_samples(std::vector<success_sample>& output)
{
    output.push_back(
        success_sample("fractional_negative", {0x66, 0x66, 0x86, 0xc0}, -4.2f, protocol_field_type::float_)
    );
    output.push_back(
        success_sample("fractional_positive", {0x66, 0x66, 0x86, 0x40}, 4.2f, protocol_field_type::float_)
    );
    output.push_back(success_sample(
        "positive_exp_positive_fractional",
        {0x01, 0x2d, 0x88, 0x61},
        3.14e20f,
        protocol_field_type::float_
    ));
    output.push_back(success_sample("zero", {0x00, 0x00, 0x00, 0x00}, 0.0f, protocol_field_type::float_));
}

void add_double_samples(std::vector<success_sample>& output)
{
    output.push_back(success_sample(
        "fractional_negative",
        {0xcd, 0xcc, 0xcc, 0xcc, 0xcc, 0xcc, 0x10, 0xc0},
        -4.2,
        protocol_field_type::double_
    ));
    output.push_back(success_sample(
        "fractional_positive",
        {0xcd, 0xcc, 0xcc, 0xcc, 0xcc, 0xcc, 0x10, 0x40},
        4.2,
        protocol_field_type::double_
    ));
    output.push_back(success_sample(
        "positive_exp_positive_fractional",
        {0xce, 0x46, 0x3c, 0x76, 0x9c, 0x68, 0x90, 0x69},
        3.14e200,
        protocol_field_type::double_
    ));
    output.push_back(success_sample(
        "zero",
        {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
        0.0,
        protocol_field_type::double_
    ));
}

void add_date_samples(std::vector<success_sample>& output)
{
    output.push_back(success_sample(
        "regular",
        {0x04, 0xda, 0x07, 0x03, 0x1c},
        date(2010u, 3u, 28u),
        protocol_field_type::date
    ));
    output.push_back(
        success_sample("min", {0x04, 0x00, 0x00, 0x01, 0x01}, date(0u, 1u, 1u), protocol_field_type::date)
    );
    output.push_back(success_sample(
        "max",
        {0x04, 0x0f, 0x27, 0x0c, 0x1f},
        date(9999u, 12u, 31u),
        protocol_field_type::date
    ));
    output.push_back(success_sample("empty", {0x00}, date(), protocol_field_type::date));
    output.push_back(success_sample("zero", {0x04, 0x00, 0x00, 0x00, 0x00}, date(), protocol_field_type::date)
    );
    output.push_back(success_sample(
        "zero_month",
        {0x04, 0xda, 0x07, 0x00, 0x01},
        date(2010u, 0u, 1u),
        protocol_field_type::date
    ));
    output.push_back(success_sample(
        "zero_day",
        {0x04, 0xda, 0x07, 0x01, 0x00},
        date(2010u, 1u, 0u),
        protocol_field_type::date
    ));
    output.push_back(success_sample(
        "zero_month_day",
        {0x04, 0xda, 0x07, 0x00, 0x00},
        date(2010u, 0u, 0u),
        protocol_field_type::date
    ));
    output.push_back(success_sample(
        "invalid_date",
        {0x04, 0xda, 0x07, 0x0b, 0x1f},
        date(2010u, 11u, 31u),
        protocol_field_type::date
    ));
}

void add_datetime_samples(protocol_field_type type, std::vector<success_sample>& output)
{
    output.push_back(
        success_sample("only_date", {0x04, 0xda, 0x07, 0x01, 0x01}, datetime(2010u, 1u, 1u), type)
    );
    output.push_back(success_sample(
        "date_h",
        {0x07, 0xda, 0x07, 0x01, 0x01, 0x14, 0x00, 0x00},
        datetime(2010u, 1u, 1u, 20u, 0u, 0u, 0u),
        type
    ));
    output.push_back(success_sample(
        "date_m",
        {0x07, 0xda, 0x07, 0x01, 0x01, 0x00, 0x01, 0x00},
        datetime(2010u, 1u, 1u, 0u, 1u, 0u, 0u),
        type
    ));
    output.push_back(success_sample(
        "date_hm",
        {0x07, 0xda, 0x07, 0x01, 0x01, 0x03, 0x02, 0x00},
        datetime(2010u, 1u, 1u, 3u, 2u, 0u, 0u),
        type
    ));
    output.push_back(success_sample(
        "date_s",
        {0x07, 0xda, 0x07, 0x01, 0x01, 0x00, 0x00, 0x01},
        datetime(2010u, 1u, 1u, 0u, 0u, 1u, 0u),
        type
    ));
    output.push_back(success_sample(
        "date_ms",
        {0x07, 0xda, 0x07, 0x01, 0x01, 0x00, 0x3b, 0x01},
        datetime(2010u, 1u, 1u, 0u, 59u, 1u, 0u),
        type
    ));
    output.push_back(success_sample(
        "date_hs",
        {0x07, 0xda, 0x07, 0x01, 0x01, 0x05, 0x00, 0x01},
        datetime(2010u, 1u, 1u, 5u, 0u, 1u, 0u),
        type
    ));
    output.push_back(success_sample(
        "date_hms",
        {0x07, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 0x3b},
        datetime(2010u, 1u, 1u, 23u, 1u, 59u, 0u),
        type
    ));
    output.push_back(success_sample(
        "date_u",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x00, 0x00, 0x00, 0x78, 0xd4, 0x03, 0x00},
        datetime(2010u, 1u, 1u, 0u, 0u, 0u, 251000u),
        type
    ));
    output.push_back(success_sample(
        "date_hu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x00, 0x00, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 1u, 1u, 23u, 0u, 0u, 967510u),
        type
    ));
    output.push_back(success_sample(
        "date_mu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x00, 0x01, 0x00, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 1u, 1u, 0u, 1u, 0u, 967510u),
        type
    ));
    output.push_back(success_sample(
        "date_hmu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 0x00, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 1u, 1u, 23u, 1u, 0u, 967510u),
        type
    ));
    output.push_back(success_sample(
        "date_su",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x00, 0x00, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 1u, 1u, 0u, 0u, 59u, 967510u),
        type
    ));
    output.push_back(success_sample(
        "date_msu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x00, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 1u, 1u, 0u, 1u, 59u, 967510u),
        type
    ));
    output.push_back(success_sample(
        "date_hsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x00, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 1u, 1u, 23u, 0u, 59u, 967510u),
        type
    ));
    output.push_back(success_sample(
        "date_hmsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 1u, 1u, 23u, 1u, 59u, 967510u),
        type
    ));

    // Invalid datetimes (because their date is invalid)
    output.push_back(success_sample("empty", {0x00}, datetime(), type));
    output.push_back(success_sample("only_date_zeros", {0x04, 0x00, 0x00, 0x00, 0x00}, datetime(), type));
    output.push_back(success_sample(
        "only_date_invalid_date",
        {0x04, 0xda, 0x07, 0x0b, 0x1f},
        datetime(2010u, 11u, 31u),
        type
    ));
    output.push_back(
        success_sample("only_date_zero_month", {0x04, 0xda, 0x07, 0x00, 0x01}, datetime(2010u, 0u, 1u), type)
    );
    output.push_back(
        success_sample("only_date_zero_day", {0x04, 0xda, 0x07, 0x01, 0x00}, datetime(2010u, 1u, 0u), type)
    );
    output.push_back(success_sample(
        "only_date_zero_month_day",
        {0x04, 0xda, 0x07, 0x00, 0x00},
        datetime(2010u, 0u, 0u),
        type
    ));

    output.push_back(
        success_sample("date_hms_zeros", {0x07, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, datetime(), type)
    );
    output.push_back(success_sample(
        "date_hms_invalid_date",
        {0x07, 0xda, 0x07, 0x0b, 0x1f, 0x17, 0x01, 0x3b},
        datetime(2010u, 11u, 31u, 23u, 1u, 59u),
        type
    ));
    output.push_back(success_sample(
        "date_hms_zero_month",
        {0x07, 0xda, 0x07, 0x00, 0x01, 0x17, 0x01, 0x3b},
        datetime(2010u, 0u, 1u, 23u, 1u, 59u),
        type
    ));
    output.push_back(success_sample(
        "date_hms_zero_day",
        {0x07, 0xda, 0x07, 0x01, 0x00, 0x17, 0x01, 0x3b},
        datetime(2010u, 1u, 0u, 23u, 1u, 59u),
        type
    ));
    output.push_back(success_sample(
        "date_hms_zero_month_day",
        {0x07, 0xda, 0x07, 0x00, 0x00, 0x17, 0x01, 0x3b},
        datetime(2010u, 0u, 0u, 23u, 1u, 59u),
        type
    ));

    output.push_back(success_sample(
        "date_hmsu_zeros",
        {0x0b, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
        datetime(),
        type
    ));
    output.push_back(success_sample(
        "date_hmsu_invalid_date",
        {0x0b, 0xda, 0x07, 0x0b, 0x1f, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 11u, 31u, 23u, 1u, 59u, 967510u),
        type
    ));
    output.push_back(success_sample(
        "date_hmsu_zero_month",
        {0x0b, 0xda, 0x07, 0x00, 0x01, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 0u, 1u, 23u, 1u, 59u, 967510u),
        type
    ));
    output.push_back(success_sample(
        "date_hmsu_zero_day",
        {0x0b, 0xda, 0x07, 0x01, 0x00, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 1u, 0u, 23u, 1u, 59u, 967510u),
        type
    ));
    output.push_back(success_sample(
        "date_hmsu_zero_month_day",
        {0x0b, 0xda, 0x07, 0x00, 0x00, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 0u, 0u, 23u, 1u, 59u, 967510u),
        type
    ));
}

void add_time_samples(std::vector<success_sample>& output)
{
    output.push_back(success_sample("zero", {0x00}, maket(0, 0, 0), protocol_field_type::time));
    output.push_back(success_sample(
        "positive_d",
        {0x08, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
        maket(48, 0, 0, 0),
        protocol_field_type::time
    ));
    output.push_back(success_sample(
        "positive_h",
        {0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x15, 0x00, 0x00},
        maket(21, 0, 0, 0),
        protocol_field_type::time
    ));
    output.push_back(success_sample(
        "positive_m",
        {0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x28, 0x00},
        maket(0, 40, 0),
        protocol_field_type::time
    ));
    output.push_back(success_sample(
        "positive_s",
        {0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x15},
        maket(0, 0, 21),
        protocol_field_type::time
    ));
    output.push_back(success_sample(
        "positive_u",
        {0x0c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xe8, 0xe5, 0x04, 0x00},
        maket(0, 0, 0, 321000),
        protocol_field_type::time
    ));
    output.push_back(success_sample(
        "positive_hmsu",
        {0x0c, 0x00, 0x22, 0x00, 0x00, 0x00, 0x16, 0x3b, 0x3a, 0x58, 0x3e, 0x0f, 0x00},
        maket(838, 59, 58, 999000),
        protocol_field_type::time
    ));
    output.push_back(success_sample(
        "negative_d",
        {0x08, 0x01, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
        -maket(48, 0, 0, 0),
        protocol_field_type::time
    ));
    output.push_back(success_sample(
        "negative_h",
        {0x08, 0x01, 0x00, 0x00, 0x00, 0x00, 0x15, 0x00, 0x00},
        -maket(21, 0, 0, 0),
        protocol_field_type::time
    ));
    output.push_back(success_sample(
        "negative_m",
        {0x08, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x28, 0x00},
        -maket(0, 40, 0),
        protocol_field_type::time
    ));
    output.push_back(success_sample(
        "negative_s",
        {0x08, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x15},
        -maket(0, 0, 21),
        protocol_field_type::time
    ));
    output.push_back(success_sample(
        "negative_u",
        {0x0c, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xe8, 0xe5, 0x04, 0x00},
        -maket(0, 0, 0, 321000),
        protocol_field_type::time
    ));
    output.push_back(success_sample(
        "negative_hmsu",
        {0x0c, 0x01, 0x22, 0x00, 0x00, 0x00, 0x16, 0x3b, 0x3a, 0x58, 0x3e, 0x0f, 0x00},
        -maket(838, 59, 58, 999000),
        protocol_field_type::time
    ));
    output.push_back(success_sample(
        "negative_sign_not_one",
        {0x0c, 0x03, 0x22, 0x00, 0x00, 0x00, 0x16, 0x3b, 0x3a, 0x58, 0x3e, 0x0f, 0x00},
        -maket(838, 59, 58, 999000),
        protocol_field_type::time
    ));
}

std::vector<success_sample> make_all_samples()
{
    std::vector<success_sample> res;
    add_string_samples(res);
    add_blob_samples(res);
    add_int_samples(res);
    add_bit_types(res);
    add_float_samples(res);
    add_double_samples(res);
    add_date_samples(res);
    add_datetime_samples(protocol_field_type::datetime, res);
    add_datetime_samples(protocol_field_type::timestamp, res);
    add_time_samples(res);
    return res;
}

BOOST_DATA_TEST_CASE(test_deserialize_binary_value_ok, data::make(make_all_samples()))
{
    auto meta = create_meta(sample.type, sample.flags, 0, sample.collation);
    field_view actual_value;
    const bytestring& buffer = sample.from;
    deserialization_context ctx(buffer.data(), buffer.data() + buffer.size(), capabilities());

    auto err = deserialize_binary_field(ctx, meta, buffer.data(), actual_value);
    BOOST_TEST(err == deserialize_errc::ok);

    // Strings are representd as string view offsets. Strings are prefixed
    // by their length, so they don't start at offset 0
    if (sample.expected.is_string() || sample.expected.is_blob())
    {
        std::size_t expected_size = sample.expected.is_string() ? sample.expected.get_string().size()
                                                                : sample.expected.get_blob().size();
        field_view expected_offset = make_svoff_fv(
            buffer.size() - expected_size,
            expected_size,
            sample.expected.is_blob()
        );
        BOOST_TEST(actual_value == expected_offset);
        field_view_access::offset_to_string_view(actual_value, buffer.data());
    }

    BOOST_TEST(actual_value == sample.expected);
    BOOST_TEST(ctx.first() == buffer.data() + buffer.size());  // all bytes consumed
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(error)

struct error_sample
{
    std::string name;
    bytestring from;
    protocol_field_type type;
    std::uint16_t flags;
    deserialize_errc expected_err;

    error_sample(
        std::string&& name,
        bytestring&& from,
        protocol_field_type type,
        std::uint16_t flags = 0,
        deserialize_errc expected_err = deserialize_errc::protocol_value_error
    )
        : name(std::move(name)), from(std::move(from)), type(type), flags(flags), expected_err(expected_err)
    {
    }

    error_sample(
        std::string&& name,
        bytestring&& from,
        protocol_field_type type,
        deserialize_errc expected_err
    )
        : name(std::move(name)), from(std::move(from)), type(type), flags(0), expected_err(expected_err)
    {
    }
};

std::ostream& operator<<(std::ostream& os, const error_sample& input)
{
    return os << "(type=" << input.type << ", name=" << input.name << ")";
}

void add_int_samples(protocol_field_type type, unsigned num_bytes, std::vector<error_sample>& output)
{
    output.emplace_back(error_sample(
        "signed_not_enough_space",
        bytestring(num_bytes, 0x0a),
        type,
        deserialize_errc::incomplete_message
    ));
    output.emplace_back(error_sample(
        "unsigned_not_enough_space",
        bytestring(num_bytes, 0x0a),
        type,
        column_flags::unsigned_,
        deserialize_errc::incomplete_message
    ));
}

void add_bit_samples(std::vector<error_sample>& output)
{
    output.emplace_back(error_sample(
        "bit_error_deserializing_string_view",
        {0x01},
        protocol_field_type::bit,
        column_flags::unsigned_,
        deserialize_errc::incomplete_message
    ));
    output.emplace_back(
        error_sample("bit_string_view_too_short", {0x00}, protocol_field_type::bit, column_flags::unsigned_)
    );
    output.emplace_back(error_sample(
        "bit_string_view_too_long",
        {0x09, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09},
        protocol_field_type::bit,
        column_flags::unsigned_
    ));
}

void add_float_samples(std::vector<error_sample>& output)
{
    output.push_back(error_sample(
        "not_enough_space",
        {0x01, 0x02, 0x03},
        protocol_field_type::float_,
        deserialize_errc::incomplete_message
    ));
    output.push_back(error_sample("inf", {0x00, 0x00, 0x80, 0x7f}, protocol_field_type::float_));
    output.push_back(error_sample("minus_inf", {0x00, 0x00, 0x80, 0xff}, protocol_field_type::float_));
    output.push_back(error_sample("nan", {0xff, 0xff, 0xff, 0x7f}, protocol_field_type::float_));
    output.push_back(error_sample("minus_nan", {0xff, 0xff, 0xff, 0xff}, protocol_field_type::float_));
}

void add_double_samples(std::vector<error_sample>& output)
{
    output.push_back(error_sample(
        "not_enough_space",
        {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07},
        protocol_field_type::double_,
        deserialize_errc::incomplete_message
    ));
    output.push_back(
        error_sample("inf", {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x7f}, protocol_field_type::double_)
    );
    output.push_back(error_sample(
        "minus_inf",
        {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xff},
        protocol_field_type::double_
    ));
    output.push_back(
        error_sample("nan", {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f}, protocol_field_type::double_)
    );
    output.push_back(error_sample(
        "minus_nan",
        {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff},
        protocol_field_type::double_
    ));
}

// Based on correct, regular date {0x04, 0xda, 0x07, 0x03, 0x1c}
void add_date_samples(std::vector<error_sample>& output)
{
    output.push_back(
        error_sample("empty", {}, protocol_field_type::date, deserialize_errc::incomplete_message)
    );
    output.push_back(error_sample(
        "incomplete_year",
        {0x04, 0xda},
        protocol_field_type::date,
        deserialize_errc::incomplete_message
    ));
    output.push_back(error_sample(
        "no_month_day",
        {0x04, 0xda, 0x07},
        protocol_field_type::date,
        deserialize_errc::incomplete_message
    ));
    output.push_back(error_sample(
        "no_day",
        {0x04, 0xda, 0x07, 0x03},
        protocol_field_type::date,
        deserialize_errc::incomplete_message
    ));
    output.push_back(error_sample(
        "invalid_year",
        {0x04, 0x10, 0x27, 0x03, 0x1c},  // year 10000
        protocol_field_type::date
    ));
    output.push_back(
        error_sample("invalid_year_max", {0x04, 0xff, 0xff, 0x03, 0x1c}, protocol_field_type::date)
    );
    output.push_back(error_sample("invalid_month", {0x04, 0xda, 0x07, 13, 0x1c}, protocol_field_type::date));
    output.push_back(
        error_sample("invalid_month_max", {0x04, 0xda, 0x07, 0xff, 0x1c}, protocol_field_type::date)
    );
    output.push_back(error_sample("invalid_day", {0x04, 0xda, 0x07, 0x03, 32}, protocol_field_type::date));
    output.push_back(
        error_sample("invalid_day_max", {0x04, 0xda, 0x07, 0x03, 0xff}, protocol_field_type::date)
    );
    output.push_back(error_sample("protocol_max", {0xff, 0xff, 0xff, 0xff, 0xff}, protocol_field_type::date));
}

// Based on correct datetime {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e,
// 0x00}
void add_datetime_samples(protocol_field_type type, std::vector<error_sample>& output)
{
    output.push_back(error_sample("empty", {}, type, deserialize_errc::incomplete_message));
    output.push_back(
        error_sample("incomplete_date", {0x04, 0xda, 0x07, 0x01}, type, deserialize_errc::incomplete_message)
    );
    output.push_back(error_sample(
        "no_hours_mins_secs",
        {0x07, 0xda, 0x07, 0x01, 0x01},
        type,
        deserialize_errc::incomplete_message
    ));
    output.push_back(error_sample(
        "no_mins_secs",
        {0x07, 0xda, 0x07, 0x01, 0x01, 0x17},
        type,
        deserialize_errc::incomplete_message
    ));
    output.push_back(error_sample(
        "no_secs",
        {0x07, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01},
        type,
        deserialize_errc::incomplete_message
    ));
    output.push_back(error_sample(
        "incomplete_micros",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e},
        type,
        deserialize_errc::incomplete_message
    ));
    output.push_back(error_sample("invalid_year_d", {0x04, 0x10, 0x27, 0x01, 0x01}, type));  // year 10000
    output.push_back(error_sample("invalid_year_hms", {0x07, 0x10, 0x27, 0x01, 0x01, 0x17, 0x01, 0x3b}, type)
    );
    output.push_back(error_sample(
        "invalid_year_hmsu",
        {0x0b, 0x10, 0x27, 0x01, 0x01, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        type
    ));
    output.push_back(error_sample(
        "invalid_year_max_hmsu",
        {0x0b, 0xff, 0xff, 0x01, 0x01, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        type
    ));
    output.push_back(error_sample("invalid_hour_hms", {0x07, 0xda, 0x07, 0x01, 0x01, 24, 0x01, 0x3b}, type));
    output.push_back(error_sample(
        "invalid_hour_hmsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 24, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        type
    ));
    output.push_back(error_sample(
        "invalid_hour_max_hmsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0xff, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        type
    ));
    output.push_back(error_sample("invalid_min_hms", {0x07, 0xda, 0x07, 0x01, 0x01, 0x17, 60, 0x3b}, type));
    output.push_back(error_sample(
        "invalid_min_hmsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 60, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        type
    ));
    output.push_back(error_sample(
        "invalid_min_max_hmsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0xff, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        type
    ));
    output.push_back(error_sample("invalid_sec_hms", {0x07, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 60}, type));
    output.push_back(error_sample(
        "invalid_sec_hmsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 60, 0x56, 0xc3, 0x0e, 0x00},
        type
    ));
    output.push_back(error_sample(
        "invalid_sec_max_hmsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 0xff, 0x56, 0xc3, 0x0e, 0x00},
        type
    ));
    output.push_back(error_sample(
        "invalid_micro_hmsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 0x3b, 0x40, 0x42, 0xf4, 0x00},
        type
    ));  // 1M
    output.push_back(error_sample(
        "invalid_micro_max_hmsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 0x3b, 0xff, 0xff, 0xff, 0xff},
        type
    ));
    output.push_back(error_sample(
        "invalid_hour_invalid_date",
        {0x0b, 0x00, 0x00, 0x00, 0x00, 0xff, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        type
    ));
    output.push_back(error_sample(
        "invalid_min_invalid_date",
        {0x0b, 0x00, 0x00, 0x00, 0x00, 0x17, 0xff, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        type
    ));
    output.push_back(error_sample(
        "invalid_sec_invalid_date",
        {0x0b, 0x00, 0x00, 0x00, 0x00, 0x17, 0x01, 0xff, 0x56, 0xc3, 0x0e, 0x00},
        type
    ));
    output.push_back(error_sample(
        "invalid_micro_invalid_date",
        {0x0b, 0x00, 0x00, 0x00, 0x00, 0x17, 0x01, 0x3b, 0xff, 0xff, 0xff, 0xff},
        type
    ));
    output.push_back(error_sample(
        "protocol_max",
        {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff},
        type
    ));
}

void add_time_samples(std::vector<error_sample>& output)
{
    constexpr auto type = protocol_field_type::time;
    output.push_back(error_sample("empty", {}, type, deserialize_errc::incomplete_message));
    output.push_back(
        error_sample("no_sign_days_hours_mins_secs", {0x08}, type, deserialize_errc::incomplete_message)
    );
    output.push_back(
        error_sample("no_days_hours_mins_secs", {0x08, 0x01}, type, deserialize_errc::incomplete_message)
    );
    output.push_back(error_sample(
        "no_hours_mins_secs",
        {0x08, 0x01, 0x22, 0x00, 0x00, 0x00},
        type,
        deserialize_errc::incomplete_message
    ));
    output.push_back(error_sample(
        "no_mins_secs",
        {0x08, 0x01, 0x22, 0x00, 0x00, 0x00, 0x16},
        type,
        deserialize_errc::incomplete_message
    ));
    output.push_back(error_sample(
        "no_secs",
        {0x08, 0x01, 0x22, 0x00, 0x00, 0x00, 0x16, 0x3b},
        type,
        deserialize_errc::incomplete_message
    ));
    output.push_back(error_sample(
        "no_micros",
        {0x0c, 0x01, 0x22, 0x00, 0x00, 0x00, 0x16, 0x3b, 0x3a},
        type,
        deserialize_errc::incomplete_message
    ));

    std::pair<const char*, std::vector<std::uint8_t>> out_of_range_cases[]{
        {"invalid_days",       {0x08, 0x00, 35, 0x00, 0x00, 0x00, 0x16, 0x3b, 0x3a}                          },
        {"invalid_days_max",   {0x08, 0x00, 0xff, 0xff, 0xff, 0xff, 0x16, 0x3b, 0x3a}                        },
        {"invalid_hours",      {0x08, 0x01, 0x22, 0x00, 0x00, 0x00, 24, 0x3b, 0x3a}                          },
        {"invalid_hours_max",  {0x08, 0x01, 0x22, 0x00, 0x00, 0x00, 0xff, 0x3b, 0x3a}                        },
        {"invalid_mins",       {0x08, 0x01, 0x22, 0x00, 0x00, 0x00, 0x16, 60, 0x3a}                          },
        {"invalid_mins_max",   {0x08, 0x01, 0x22, 0x00, 0x00, 0x00, 0x16, 0xff, 0x3a}                        },
        {"invalid_secs",       {0x08, 0x01, 0x22, 0x00, 0x00, 0x00, 0x16, 0x3b, 60}                          },
        {"invalid_secs_max",   {0x08, 0x01, 0x22, 0x00, 0x00, 0x00, 0x16, 0x3b, 0xff}                        },
        {"invalid_micros",     {0x0c, 0x01, 0x22, 0x00, 0x00, 0x00, 0x16, 0x3b, 0x3a, 0x40, 0x42, 0xf4, 0x00}},
        {"invalid_micros_max",
         {0x0c, 0x01, 0x22, 0x00, 0x00, 0x00, 0x16, 0x3b, 0x3a, 0xff, 0xff, 0xff, 0xff}                      },
    };

    for (auto& c : out_of_range_cases)
    {
        // Positive
        c.second[1] = 0x00;
        output.emplace_back(c.first + std::string("_positive"), bytestring(c.second), type);

        // Negative
        c.second[1] = 0x01;
        output.emplace_back(c.first + std::string("_negative"), std::move(c.second), type);
    }
}

std::vector<error_sample> make_all_samples()
{
    std::vector<error_sample> res;
    add_int_samples(protocol_field_type::tiny, 0, res);
    add_int_samples(protocol_field_type::short_, 1, res);
    add_int_samples(protocol_field_type::int24, 3, res);
    add_int_samples(protocol_field_type::long_, 3, res);
    add_int_samples(protocol_field_type::longlong, 7, res);
    add_int_samples(protocol_field_type::year, 1, res);
    add_bit_samples(res);
    add_float_samples(res);
    add_double_samples(res);
    add_date_samples(res);
    add_datetime_samples(protocol_field_type::datetime, res);
    add_datetime_samples(protocol_field_type::timestamp, res);
    add_time_samples(res);
    return res;
}

BOOST_DATA_TEST_CASE(test_deserialize_binary_value_error, data::make(make_all_samples()))
{
    auto meta = create_meta(sample.type, sample.flags);
    field_view actual_value;
    const bytestring& buff = sample.from;
    deserialization_context ctx(buff.data(), buff.data() + buff.size(), capabilities());
    auto err = deserialize_binary_field(ctx, meta, buff.data(), actual_value);
    auto expected = sample.expected_err;
    BOOST_TEST(expected == err);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

}  // namespace
