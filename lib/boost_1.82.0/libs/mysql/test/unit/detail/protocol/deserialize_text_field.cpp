//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

// Test deserialize_text_value(), just positive cases

#include <boost/mysql/blob_view.hpp>
#include <boost/mysql/date.hpp>
#include <boost/mysql/datetime.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/mysql_collations.hpp>

#include <boost/mysql/detail/auxiliar/access_fwd.hpp>
#include <boost/mysql/detail/auxiliar/string_view_offset.hpp>
#include <boost/mysql/detail/auxiliar/stringize.hpp>
#include <boost/mysql/detail/protocol/constants.hpp>
#include <boost/mysql/detail/protocol/deserialize_text_field.hpp>

#include <boost/test/data/monomorphic/collection.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test_suite.hpp>

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
using boost::mysql::string_view;

namespace {

// clang-format off

BOOST_AUTO_TEST_SUITE(test_deserialize_text_field)

// Success cases
BOOST_AUTO_TEST_SUITE(success)

struct success_sample
{
    std::string name;
    std::string from;
    field_view expected;
    protocol_field_type type;
    std::uint16_t flags;
    unsigned decimals;
    std::uint16_t collation;

    template <class T>
    success_sample(
        std::string name,
        std::string from,
        T&& expected_value,
        protocol_field_type type,
        std::uint16_t flags=0,
        unsigned decimals=0,
        std::uint16_t collation=boost::mysql::mysql_collations::utf8mb4_general_ci 
    ) :
        name(std::move(name)),
        from(std::move(from)),
        expected(std::forward<T>(expected_value)),
        type(type),
        flags(flags),
        decimals(decimals),
        collation(collation)
    {
    }
};

std::ostream& operator<<(std::ostream& os, const success_sample& input)
{
    return os << "(input=" << input.from
              << ", type=" << input.type
              << ", name=" << input.name
              << ")";
}

void add_string_samples(std::vector<success_sample>& output)
{
    output.emplace_back("varchar_non_empty", "string", "string", protocol_field_type::var_string);
    output.emplace_back("varchar_empty", "", "", protocol_field_type::var_string);
    output.emplace_back("char", "", "", protocol_field_type::string);
    output.emplace_back("text", "value", "value", protocol_field_type::blob);
    output.emplace_back("enum", "value", "value", protocol_field_type::string, column_flags::enum_);
    output.emplace_back("set", "value1,value2", "value1,value2", protocol_field_type::string, column_flags::set);
    output.emplace_back("decimal", "1", "1", protocol_field_type::newdecimal);
    output.emplace_back("json", "{}", "{}", protocol_field_type::json);
}

void add_blob_samples(std::vector<success_sample>& output)
{
    static constexpr std::uint8_t buff [] = { 0x00, 0x01, 0x02, 0x03 };
    std::string from { 0x00, 0x01, 0x02, 0x03 };

    output.emplace_back("varbinary_non_empty", from, blob_view(buff), protocol_field_type::var_string, column_flags::binary, 0, binary_collation);
    output.emplace_back("varbinary_empty", "", blob_view(), protocol_field_type::var_string, column_flags::binary, 0, binary_collation);
    output.emplace_back("binary", from, blob_view(buff), protocol_field_type::string, column_flags::binary, 0, binary_collation);
    output.emplace_back("blob", from, blob_view(buff), protocol_field_type::blob, column_flags::binary, 0, binary_collation);
    output.emplace_back("geometry", from, blob_view(buff), protocol_field_type::geometry,
            column_flags::binary | column_flags::blob, 0, binary_collation);

    // Anything we don't know what it is, we interpret as a blob
    output.emplace_back("unknown_protocol_type", from,
            blob_view(buff), static_cast<protocol_field_type>(0x23));
}

void add_int_samples_helper(
    std::string signed_max_s,
    std::int64_t signed_max_b,
    std::string signed_min_s,
    std::int64_t signed_min_b,
    std::string unsigned_max_s,
    std::uint64_t unsigned_max_b,
    std::string zerofill_s,
    std::uint64_t zerofill_b,
    protocol_field_type type,
    std::vector<success_sample>& output
)
{
    output.emplace_back("signed", "20", std::int64_t(20), type);
    output.emplace_back("signed_max", std::move(signed_max_s), signed_max_b, type);
    output.emplace_back("signed_negative", "-20", std::int64_t(-20), type);
    output.emplace_back("signed_min", std::move(signed_min_s), signed_min_b, type);
    output.emplace_back("unsigned", "20", std::uint64_t(20), type, column_flags::unsigned_);
    output.emplace_back("unsigned_min", "0", std::uint64_t(0), type, column_flags::unsigned_);
    output.emplace_back("unsigned_max", std::move(unsigned_max_s),
        unsigned_max_b, type, column_flags::unsigned_);
    output.emplace_back("unsigned_zerofill", std::move(zerofill_s),
        zerofill_b, type, column_flags::unsigned_ | column_flags::zerofill);
}

void add_int_samples(std::vector<success_sample>& output)
{
    add_int_samples_helper(
        "127", 127,
        "-128", -128,
        "255", 255,
        "010", 10,
        protocol_field_type::tiny,
        output
    );
    add_int_samples_helper(
        "32767", 32767,
        "-32768", -32768,
        "65535", 65535,
        "00535", 535,
        protocol_field_type::short_,
        output
    );
    add_int_samples_helper(
        "8388607", 8388607,
        "-8388608",-8388608,
        "16777215", 16777215,
        "00007215", 7215,
        protocol_field_type::int24,
        output
    );
    add_int_samples_helper(
        "2147483647", 2147483647,
        "-2147483648", -2147483647 - 1, // minus is not part of literal, avoids warning
        "4294967295", 4294967295,
        "0000067295", 67295,
        protocol_field_type::long_,
        output
    );
    add_int_samples_helper(
        "9223372036854775807", 9223372036854775807,
        "-9223372036854775808", -9223372036854775807 - 1, // minus is not part of literal, avoids warning
        "18446744073709551615", 18446744073709551615ULL, // suffix needed to avoid warning
        "0000067295", 67295,
        protocol_field_type::longlong,
        output
    );
    add_int_samples_helper(
        "9223372036854775807", 9223372036854775807,
        "-9223372036854775808", -9223372036854775807 - 1, // minus is not part of literal, avoids warning
        "18446744073709551615", 18446744073709551615ULL, // suffix needed to avoid warning
        "0000067295", 67295,
        protocol_field_type::longlong,
        output
    );

    // YEAR
    output.emplace_back("regular_value", "1999", std::uint64_t(1999), protocol_field_type::year, column_flags::unsigned_);
    output.emplace_back("min", "1901", std::uint64_t(1901), protocol_field_type::year, column_flags::unsigned_);
    output.emplace_back("max", "2155", std::uint64_t(2155), protocol_field_type::year, column_flags::unsigned_);
    output.emplace_back("zero", "0000", std::uint64_t(0), protocol_field_type::year, column_flags::unsigned_);
}

// bit
void add_bit_types(std::vector<success_sample>& output)
{
    output.emplace_back("bit_8", "\x12", std::uint64_t(0x12), 
        protocol_field_type::bit, column_flags::unsigned_);
    output.emplace_back("bit_16", "\x12\x34", std::uint64_t(0x1234), 
        protocol_field_type::bit, column_flags::unsigned_);
    output.emplace_back("bit_24", "\x12\x34\x56", std::uint64_t(0x123456), 
        protocol_field_type::bit, column_flags::unsigned_);
    output.emplace_back("bit_32", "\x12\x34\x56\x78", std::uint64_t(0x12345678), 
        protocol_field_type::bit, column_flags::unsigned_);
    output.emplace_back("bit_40", "\x12\x34\x56\x78\x9a", std::uint64_t(0x123456789a), 
        protocol_field_type::bit, column_flags::unsigned_);
    output.emplace_back("bit_48", "\x12\x34\x56\x78\x9a\xbc", std::uint64_t(0x123456789abc), 
        protocol_field_type::bit, column_flags::unsigned_);
    output.emplace_back("bit_56", "\x12\x34\x56\x78\x9a\xbc\xde", std::uint64_t(0x123456789abcde), 
        protocol_field_type::bit, column_flags::unsigned_);
    output.emplace_back("bit_64", "\x12\x34\x56\x78\x9a\xbc\xde\xf0", std::uint64_t(0x123456789abcdef0), 
        protocol_field_type::bit, column_flags::unsigned_);
}

template <class T>
void add_float_samples(
    protocol_field_type type,
    std::vector<success_sample>& output
)
{
    output.emplace_back("zero", "0", T(0.0),  type);
    output.emplace_back("integer_positive", "4", T(4.0),  type);
    output.emplace_back("integer_negative", "-5", T(-5.0),  type);
    output.emplace_back("fractional_positive", "3.147", T(3.147),  type);
    output.emplace_back("fractional_negative", "-3.147", T(-3.147),  type);
    output.emplace_back("positive_exponent_positive_integer", "3e20", T(3e20),  type);
    output.emplace_back("positive_exponent_negative_integer", "-3e20", T(-3e20),  type);
    output.emplace_back("positive_exponent_positive_fractional", "3.14e20", T(3.14e20),  type);
    output.emplace_back("positive_exponent_negative_fractional", "-3.45e20", T(-3.45e20),  type);
    output.emplace_back("negative_exponent_positive_integer", "3e-20", T(3e-20),  type);
    output.emplace_back("negative_exponent_negative_integer", "-3e-20", T(-3e-20),  type);
    output.emplace_back("negative_exponent_positive_fractional", "3.14e-20", T(3.14e-20),  type);
    output.emplace_back("negative_exponent_negative_fractional", "-3.45e-20", T(-3.45e-20),  type);
}

void add_date_samples(std::vector<success_sample>& output)
{
    output.emplace_back("regular_date", "2019-02-28", date(2019u, 2u, 28u), protocol_field_type::date);
    output.emplace_back("leap_year", "1788-02-29", date(1788u, 2u, 29u), protocol_field_type::date);
    output.emplace_back("min", "0000-01-01", date(0u, 1u, 1u), protocol_field_type::date);
    output.emplace_back("max", "9999-12-31", date(9999u, 12u, 31u), protocol_field_type::date);
    output.emplace_back("zero", "0000-00-00", date(), protocol_field_type::date);
    output.emplace_back("zero_month", "0000-00-01", date(0u, 0u, 1u), protocol_field_type::date);
    output.emplace_back("zero_day", "0000-01-00", date(0u, 1u, 0u), protocol_field_type::date);
    output.emplace_back("zero_month_day_nonzero_year", "2010-00-00", date(2010u, 0u, 0u), protocol_field_type::date);
    output.emplace_back("invalid_date", "2010-11-31", date(2010u, 11u, 31u), protocol_field_type::date);
}

void add_datetime_samples(
    protocol_field_type type,
    std::vector<success_sample>& output
)
{
    output.emplace_back("0_decimals_date", "2010-02-15 00:00:00", datetime(2010u, 2u, 15u), type);
    output.emplace_back("0_decimals_h", "2010-02-15 02:00:00", datetime(2010u, 2u, 15u, 2u), type);
    output.emplace_back("0_decimals_hm", "2010-02-15 02:05:00", datetime(2010u, 2u, 15u, 2u, 5u), type);
    output.emplace_back("0_decimals_hms", "2010-02-15 02:05:30", datetime(2010u, 2u, 15u, 2u, 5u, 30u), type);
    output.emplace_back("0_decimals_min", "0000-01-01 00:00:00", datetime(0u, 1u, 1u), type);
    output.emplace_back("0_decimals_max", "9999-12-31 23:59:59", datetime(9999u, 12u, 31u, 23u, 59u, 59u), type);

    output.emplace_back("1_decimals_date", "2010-02-15 00:00:00.0", datetime(2010u, 2u, 15u), type, 0, 1);
    output.emplace_back("1_decimals_h", "2010-02-15 02:00:00.0", datetime(2010u, 2u, 15u, 2u), type, 0, 1);
    output.emplace_back("1_decimals_hm", "2010-02-15 02:05:00.0", datetime(2010u, 2u, 15u, 2u, 5u), type, 0, 1);
    output.emplace_back("1_decimals_hms", "2010-02-15 02:05:30.0", datetime(2010u, 2u, 15u, 2u, 5u, 30u), type, 0, 1);
    output.emplace_back("1_decimals_hmsu", "2010-02-15 02:05:30.5", datetime(2010u, 2u, 15u, 2u, 5u, 30u, 500000u), type, 0, 1);
    output.emplace_back("1_decimals_min", "0000-01-01 00:00:00.0", datetime(0u, 1u, 1u), type, 0, 1);
    output.emplace_back("1_decimals_max", "9999-12-31 23:59:59.9", datetime(9999u, 12u, 31u, 23u, 59u, 59u, 900000u), type, 0, 1);

    output.emplace_back("2_decimals_hms", "2010-02-15 02:05:30.00", datetime(2010u, 2u, 15u, 2u, 5u, 30u), type, 0, 2);
    output.emplace_back("2_decimals_hmsu", "2010-02-15 02:05:30.05", datetime(2010u, 2u, 15u, 2u, 5u, 30u, 50000u), type, 0, 2);
    output.emplace_back("2_decimals_min", "0000-01-01 00:00:00.00", datetime(0u, 1u, 1u), type, 0, 2);
    output.emplace_back("2_decimals_max", "9999-12-31 23:59:59.99", datetime(9999u, 12u, 31u, 23u, 59u, 59u, 990000u), type, 0, 2);

    output.emplace_back("3_decimals_hms", "2010-02-15 02:05:30.000", datetime(2010u, 2u, 15u, 2u, 5u, 30u), type, 0, 3);
    output.emplace_back("3_decimals_hmsu", "2010-02-15 02:05:30.420", datetime(2010u, 2u, 15u, 2u, 5u, 30u, 420000u), type, 0, 3);
    output.emplace_back("3_decimals_min", "0000-01-01 00:00:00.000", datetime(0u, 1u, 1u), type, 0, 3);
    output.emplace_back("3_decimals_max", "9999-12-31 23:59:59.999", datetime(9999u, 12u, 31u, 23u, 59u, 59u, 999000u), type, 0, 3);

    output.emplace_back("4_decimals_hms", "2010-02-15 02:05:30.0000", datetime(2010u, 2u, 15u, 2u, 5u, 30u), type, 0, 4);
    output.emplace_back("4_decimals_hmsu", "2010-02-15 02:05:30.4267", datetime(2010u, 2u, 15u, 2u, 5u, 30u, 426700u), type, 0, 4);
    output.emplace_back("4_decimals_min", "0000-01-01 00:00:00.0000", datetime(0u, 1u, 1u), type, 0, 4);
    output.emplace_back("4_decimals_max", "9999-12-31 23:59:59.9999", datetime(9999u, 12u, 31u, 23u, 59u, 59u, 999900u), type, 0, 4);

    output.emplace_back("5_decimals_hms", "2010-02-15 02:05:30.00000", datetime(2010u, 2u, 15u, 2u, 5u, 30u), type, 0, 5);
    output.emplace_back("5_decimals_hmsu", "2010-02-15 02:05:30.00239", datetime(2010u, 2u, 15u, 2u, 5u, 30u, 2390u), type, 0, 5);
    output.emplace_back("5_decimals_min", "0000-01-01 00:00:00.00000", datetime(0u, 1u, 1u), type, 0, 5);
    output.emplace_back("5_decimals_max", "9999-12-31 23:59:59.99999", datetime(9999u, 12u, 31u, 23u, 59u, 59u, 999990u), type, 0, 5);

    output.emplace_back("6_decimals_hms", "2010-02-15 02:05:30.000000", datetime(2010u, 2u, 15u, 2u, 5u, 30u), type, 0, 6);
    output.emplace_back("6_decimals_hmsu", "2010-02-15 02:05:30.002395", datetime(2010u, 2u, 15u, 2u, 5u, 30u, 2395u), type, 0, 6);
    output.emplace_back("6_decimals_min", "0000-01-01 00:00:00.000000", datetime(0u, 1u, 1u), type, 0, 6);
    output.emplace_back("6_decimals_max", "9999-12-31 23:59:59.999999", datetime(9999u, 12u, 31u, 23u, 59u, 59u, 999999u), type, 0, 6);

    // not a real case, we cap decimals to 6
    output.emplace_back("7_decimals", "2010-02-15 02:05:30.002395", datetime(2010, 2, 15, 2, 5, 30, 2395), type, 0, 7);

    // Invalid datetimes (because their date is invalid)
    output.emplace_back("0_decimals_zero", "0000-00-00 00:00:00", datetime(), type);
    output.emplace_back("0_decimals_invalid_date", "2010-11-31 01:10:59", datetime(2010u, 11u, 31u, 1u, 10u, 59u), type);
    output.emplace_back("0_decimals_zero_month", "2010-00-31 01:10:59", datetime(2010u, 0u, 31u, 1u, 10u, 59u), type);
    output.emplace_back("0_decimals_zero_day", "2010-11-00 01:10:59", datetime(2010u, 11u, 0u, 1u, 10u, 59u), type);
    output.emplace_back("0_decimals_zero_month_day", "2010-00-00 01:10:59", datetime(2010u, 0u, 0u, 1u, 10u, 59u), type);

    output.emplace_back("1_decimals_zero", "0000-00-00 00:00:00.0", datetime(), type, 0, 1);
    output.emplace_back("1_decimals_invalid_date", "2010-11-31 01:10:59.9", datetime(2010u, 11u, 31u, 1u, 10u, 59u, 900000u), type, 0, 1);
    output.emplace_back("1_decimals_zero_month", "2010-00-31 01:10:59.9", datetime(2010u, 0u, 31u, 1u, 10u, 59u, 900000u), type, 0, 1);
    output.emplace_back("1_decimals_zero_day", "2010-11-00 01:10:59.9", datetime(2010u, 11u, 0u, 1u, 10u, 59u, 900000u), type, 0, 1);
    output.emplace_back("1_decimals_zero_month_day", "2010-00-00 01:10:59.9", datetime(2010u, 0u, 0u, 1u, 10u, 59u, 900000u), type, 0, 1);

    output.emplace_back("2_decimals_zero", "0000-00-00 00:00:00.00", datetime(), type, 0, 2);
    output.emplace_back("2_decimals_invalid_date", "2010-11-31 01:10:59.98", datetime(2010u, 11u, 31u, 1u, 10u, 59u, 980000u), type, 0, 2);
    output.emplace_back("2_decimals_zero_month", "2010-00-31 01:10:59.98", datetime(2010u, 0u, 31u, 1u, 10u, 59u, 980000u), type, 0, 2);
    output.emplace_back("2_decimals_zero_day", "2010-11-00 01:10:59.98", datetime(2010u, 11u, 0u, 1u, 10u, 59u, 980000u), type, 0, 2);
    output.emplace_back("2_decimals_zero_month_day", "2010-00-00 01:10:59.98", datetime(2010u, 0u, 0u, 1u, 10u, 59u, 980000u), type, 0, 2);

    output.emplace_back("3_decimals_zero", "0000-00-00 00:00:00.000", datetime(), type, 0, 3);
    output.emplace_back("3_decimals_invalid_date", "2010-11-31 01:10:59.987", datetime(2010u, 11u, 31u, 1u, 10u, 59u, 987000u), type, 0, 3);
    output.emplace_back("3_decimals_zero_month", "2010-00-31 01:10:59.987", datetime(2010u, 0u, 31u, 1u, 10u, 59u, 987000u), type, 0, 3);
    output.emplace_back("3_decimals_zero_day", "2010-11-00 01:10:59.987", datetime(2010u, 11u, 0u, 1u, 10u, 59u, 987000u), type, 0, 3);
    output.emplace_back("3_decimals_zero_month_day", "2010-00-00 01:10:59.987", datetime(2010u, 0u, 0u, 1u, 10u, 59u, 987000u), type, 0, 3);

    output.emplace_back("4_decimals_zero", "0000-00-00 00:00:00.0000", datetime(), type, 0, 4);
    output.emplace_back("4_decimals_invalid_date", "2010-11-31 01:10:59.9876", datetime(2010u, 11u, 31u, 1u, 10u, 59u, 987600u), type, 0, 4);
    output.emplace_back("4_decimals_zero_month", "2010-00-31 01:10:59.9876", datetime(2010u, 0u, 31u, 1u, 10u, 59u, 987600u), type, 0, 4);
    output.emplace_back("4_decimals_zero_day", "2010-11-00 01:10:59.9876", datetime(2010u, 11u, 0u, 1u, 10u, 59u, 987600u), type, 0, 4);
    output.emplace_back("4_decimals_zero_month_day", "2010-00-00 01:10:59.9876", datetime(2010u, 0u, 0u, 1u, 10u, 59u, 987600u), type, 0, 4);

    output.emplace_back("5_decimals_zero", "0000-00-00 00:00:00.00000", datetime(), type, 0, 5);
    output.emplace_back("5_decimals_invalid_date", "2010-11-31 01:10:59.98765", datetime(2010u, 11u, 31u, 1u, 10u, 59u, 987650u), type, 0, 5);
    output.emplace_back("5_decimals_zero_month", "2010-00-31 01:10:59.98765", datetime(2010u, 0u, 31u, 1u, 10u, 59u, 987650u), type, 0, 5);
    output.emplace_back("5_decimals_zero_day", "2010-11-00 01:10:59.98765", datetime(2010u, 11u, 0u, 1u, 10u, 59u, 987650u), type, 0, 5);
    output.emplace_back("5_decimals_zero_month_day", "2010-00-00 01:10:59.98765", datetime(2010u, 0u, 0u, 1u, 10u, 59u, 987650u), type, 0, 5);

    output.emplace_back("6_decimals_zero", "0000-00-00 00:00:00.000000", datetime(), type, 0, 6);
    output.emplace_back("6_decimals_invalid_date", "2010-11-31 01:10:59.987654", datetime(2010u, 11u, 31u, 1u, 10u, 59u, 987654u), type, 0, 6);
    output.emplace_back("6_decimals_zero_month", "2010-00-31 01:10:59.987654", datetime(2010u, 0u, 31u, 1u, 10u, 59u, 987654u), type, 0, 6);
    output.emplace_back("6_decimals_zero_day", "2010-11-00 01:10:59.987654", datetime(2010u, 11u, 0u, 1u, 10u, 59u, 987654u), type, 0, 6);
    output.emplace_back("6_decimals_zero_month_day", "2010-00-00 01:10:59.987654", datetime(2010u, 0u, 0u, 1u, 10u, 59u, 987654u), type, 0, 6);
}

void add_time_samples(std::vector<success_sample>& output)
{
    output.emplace_back("0_decimals_positive_h", "01:00:00", maket(1, 0, 0), protocol_field_type::time);
    output.emplace_back("0_decimals_positive_hm", "12:03:00", maket(12, 3, 0), protocol_field_type::time);
    output.emplace_back("0_decimals_positive_hms", "14:51:23", maket(14, 51, 23), protocol_field_type::time);
    output.emplace_back("0_decimals_max", "838:59:59", maket(838, 59, 59), protocol_field_type::time);
    output.emplace_back("0_decimals_negative_h", "-06:00:00", -maket(6, 0, 0), protocol_field_type::time);
    output.emplace_back("0_decimals_negative_hm", "-12:03:00", -maket(12, 3, 0), protocol_field_type::time);
    output.emplace_back("0_decimals_negative_hms", "-14:51:23", -maket(14, 51, 23), protocol_field_type::time);
    output.emplace_back("0_decimals_min", "-838:59:59", -maket(838, 59, 59), protocol_field_type::time);
    output.emplace_back("0_decimals_zero", "00:00:00", maket(0, 0, 0), protocol_field_type::time);
    output.emplace_back("0_decimals_negative_h0", "-00:51:23", -maket(0, 51, 23), protocol_field_type::time);

    output.emplace_back("1_decimals_positive_hms", "14:51:23.0", maket(14, 51, 23), protocol_field_type::time, 0, 1);
    output.emplace_back("1_decimals_positive_hmsu", "14:51:23.5", maket(14, 51, 23, 500000), protocol_field_type::time, 0, 1);
    output.emplace_back("1_decimals_max", "838:59:58.9", maket(838, 59, 58, 900000), protocol_field_type::time, 0, 1);
    output.emplace_back("1_decimals_negative_hms", "-14:51:23.0", -maket(14, 51, 23), protocol_field_type::time, 0, 1);
    output.emplace_back("1_decimals_negative_hmsu", "-14:51:23.5", -maket(14, 51, 23, 500000), protocol_field_type::time, 0, 1);
    output.emplace_back("1_decimals_min", "-838:59:58.9", -maket(838, 59, 58, 900000), protocol_field_type::time, 0, 1);
    output.emplace_back("1_decimals_zero", "00:00:00.0", maket(0, 0, 0), protocol_field_type::time, 0, 1);
    output.emplace_back("1_decimals_negative_h0", "-00:51:23.1", -maket(0, 51, 23, 100000), protocol_field_type::time, 0, 1);

    output.emplace_back("2_decimals_positive_hms", "14:51:23.00", maket(14, 51, 23), protocol_field_type::time, 0, 2);
    output.emplace_back("2_decimals_positive_hmsu", "14:51:23.52", maket(14, 51, 23, 520000), protocol_field_type::time, 0, 2);
    output.emplace_back("2_decimals_max", "838:59:58.99", maket(838, 59, 58, 990000), protocol_field_type::time, 0, 2);
    output.emplace_back("2_decimals_negative_hms", "-14:51:23.00", -maket(14, 51, 23), protocol_field_type::time, 0, 2);
    output.emplace_back("2_decimals_negative_hmsu", "-14:51:23.50", -maket(14, 51, 23, 500000), protocol_field_type::time, 0, 2);
    output.emplace_back("2_decimals_min", "-838:59:58.99", -maket(838, 59, 58, 990000), protocol_field_type::time, 0, 2);
    output.emplace_back("2_decimals_zero", "00:00:00.00", maket(0, 0, 0), protocol_field_type::time, 0, 2);
    output.emplace_back("2_decimals_negative_h0", "-00:51:23.12", -maket(0, 51, 23, 120000), protocol_field_type::time, 0, 2);

    output.emplace_back("3_decimals_positive_hms", "14:51:23.000", maket(14, 51, 23), protocol_field_type::time, 0, 3);
    output.emplace_back("3_decimals_positive_hmsu", "14:51:23.501", maket(14, 51, 23, 501000), protocol_field_type::time, 0, 3);
    output.emplace_back("3_decimals_max", "838:59:58.999", maket(838, 59, 58, 999000), protocol_field_type::time, 0, 3);
    output.emplace_back("3_decimals_negative_hms", "-14:51:23.000", -maket(14, 51, 23), protocol_field_type::time, 0, 3);
    output.emplace_back("3_decimals_negative_hmsu", "-14:51:23.003", -maket(14, 51, 23, 3000), protocol_field_type::time, 0, 3);
    output.emplace_back("3_decimals_min", "-838:59:58.999", -maket(838, 59, 58, 999000), protocol_field_type::time, 0, 3);
    output.emplace_back("3_decimals_zero", "00:00:00.000", maket(0, 0, 0), protocol_field_type::time, 0, 3);
    output.emplace_back("3_decimals_negative_h0", "-00:51:23.123", -maket(0, 51, 23, 123000), protocol_field_type::time, 0, 3);

    output.emplace_back("4_decimals_positive_hms", "14:51:23.0000", maket(14, 51, 23), protocol_field_type::time, 0, 4);
    output.emplace_back("4_decimals_positive_hmsu", "14:51:23.5017", maket(14, 51, 23, 501700), protocol_field_type::time, 0, 4);
    output.emplace_back("4_decimals_max", "838:59:58.9999", maket(838, 59, 58, 999900), protocol_field_type::time, 0, 4);
    output.emplace_back("4_decimals_negative_hms", "-14:51:23.0000", -maket(14, 51, 23), protocol_field_type::time, 0, 4);
    output.emplace_back("4_decimals_negative_hmsu", "-14:51:23.0038", -maket(14, 51, 23, 3800), protocol_field_type::time, 0, 4);
    output.emplace_back("4_decimals_min", "-838:59:58.9999", -maket(838, 59, 58, 999900), protocol_field_type::time, 0, 4);
    output.emplace_back("4_decimals_zero", "00:00:00.0000", maket(0, 0, 0), protocol_field_type::time, 0, 4);
    output.emplace_back("4_decimals_negative_h0", "-00:51:23.1234", -maket(0, 51, 23, 123400), protocol_field_type::time, 0, 4);

    output.emplace_back("5_decimals_positive_hms", "14:51:23.00000", maket(14, 51, 23), protocol_field_type::time, 0, 5);
    output.emplace_back("5_decimals_positive_hmsu", "14:51:23.50171", maket(14, 51, 23, 501710), protocol_field_type::time, 0, 5);
    output.emplace_back("5_decimals_max", "838:59:58.99999", maket(838, 59, 58, 999990), protocol_field_type::time, 0, 5);
    output.emplace_back("5_decimals_negative_hms", "-14:51:23.00000", -maket(14, 51, 23), protocol_field_type::time, 0, 5);
    output.emplace_back("5_decimals_negative_hmsu", "-14:51:23.00009", -maket(14, 51, 23, 90), protocol_field_type::time, 0, 5);
    output.emplace_back("5_decimals_min", "-838:59:58.99999", -maket(838, 59, 58, 999990), protocol_field_type::time, 0, 5);
    output.emplace_back("5_decimals_zero", "00:00:00.00000", maket(0, 0, 0), protocol_field_type::time, 0, 5);
    output.emplace_back("5_decimals_negative_h0", "-00:51:23.12345", -maket(0, 51, 23, 123450), protocol_field_type::time, 0, 5);

    output.emplace_back("6_decimals_positive_hms", "14:51:23.000000", maket(14, 51, 23), protocol_field_type::time, 0, 6);
    output.emplace_back("6_decimals_positive_hmsu", "14:51:23.501717", maket(14, 51, 23, 501717), protocol_field_type::time, 0, 6);
    output.emplace_back("6_decimals_max", "838:59:58.999999", maket(838, 59, 58, 999999), protocol_field_type::time, 0, 6);
    output.emplace_back("6_decimals_negative_hms", "-14:51:23.000000", -maket(14, 51, 23), protocol_field_type::time, 0, 6);
    output.emplace_back("6_decimals_negative_hmsu", "-14:51:23.900000", -maket(14, 51, 23, 900000), protocol_field_type::time, 0, 6);
    output.emplace_back("6_decimals_min", "-838:59:58.999999", -maket(838, 59, 58, 999999), protocol_field_type::time, 0, 6);
    output.emplace_back("6_decimals_zero", "00:00:00.000000", maket(0, 0, 0), protocol_field_type::time, 0, 6);
    output.emplace_back("6_decimals_negative_h0", "-00:51:23.123456", -maket(0, 51, 23, 123456), protocol_field_type::time, 0, 6);

    // This is not a real case - we cap anything above 6 decimals to 6
    output.emplace_back("7_decimals", "14:51:23.501717", maket(14, 51, 23, 501717), protocol_field_type::time, 0, 7);
}

std::vector<success_sample> make_all_samples()
{
    std::vector<success_sample> res;
    add_string_samples(res);
    add_blob_samples(res);
    add_int_samples(res);
    add_bit_types(res);
    add_float_samples<float>(protocol_field_type::float_, res);
    add_float_samples<double>(protocol_field_type::double_, res);
    add_date_samples(res);
    add_datetime_samples(protocol_field_type::datetime, res);
    add_datetime_samples(protocol_field_type::timestamp, res);
    add_time_samples(res);
    return res;
}

BOOST_DATA_TEST_CASE(ok, data::make(make_all_samples()))
{
    auto meta = create_meta(sample.type, sample.flags, static_cast<std::uint8_t>(sample.decimals), sample.collation);
    const std::uint8_t* buffer_first = reinterpret_cast<const std::uint8_t*>(sample.from.data());
    field_view actual_value;

    auto err = deserialize_text_field(
        sample.from,
        meta,
        buffer_first,
        actual_value
    );
    BOOST_TEST(err == deserialize_errc::ok);

    // Strings are representd as string view offsets
    if (sample.expected.is_string() || sample.expected.is_blob())
    {
        std::size_t expected_offset = sample.expected.is_string() ? sample.expected.get_string().size() : sample.expected.get_blob().size();
        field_view expected_offset_fv = make_svoff_fv(0, expected_offset, sample.expected.is_blob());
        BOOST_TEST(actual_value == expected_offset_fv);
        field_view_access::offset_to_string_view(actual_value, buffer_first);
    }

    BOOST_TEST(actual_value == sample.expected);
}

BOOST_AUTO_TEST_SUITE_END()

//
// Error cases
//

BOOST_AUTO_TEST_SUITE(errors)

struct error_sample
{
    std::string name;
    string_view from;
    protocol_field_type type;
    std::uint16_t flags;
    unsigned decimals;
    deserialize_errc expected_err;

    error_sample(std::string&& name, string_view from, protocol_field_type type,
            std::uint16_t flags=0, unsigned decimals=0, deserialize_errc expected_err = deserialize_errc::protocol_value_error) :
        name(std::move(name)),
        from(from),
        type(type),
        flags(flags),
        decimals(decimals),
        expected_err(expected_err)
    {
    }
};

std::ostream& operator<<(std::ostream& os, const error_sample& input)
{
    return os << "(input=" << input.from
              << ", type=" << input.type
              << ", name=" << input.name
              << ")";
}

void add_int_samples(
    protocol_field_type t,
    std::vector<error_sample>& output
)
{
    output.emplace_back("signed_blank", "", t);
    output.emplace_back("signed_non_number", "abtrf", t);
    output.emplace_back("signed_hex", "0x01", t);
    output.emplace_back("signed_fractional", "1.1", t);
    output.emplace_back("signed_exp", "2e10", t);
    output.emplace_back("signed_lt_min", "-9223372036854775809", t);
    output.emplace_back("signed_gt_max", "9223372036854775808", t);
    output.emplace_back("unsigned_blank", "", t, column_flags::unsigned_);
    output.emplace_back("unsigned_non_number", "abtrf", t, column_flags::unsigned_);
    output.emplace_back("unsigned_hex", "0x01", t, column_flags::unsigned_);
    output.emplace_back("unsigned_fractional", "1.1", t, column_flags::unsigned_);
    output.emplace_back("unsigned_exp", "2e10", t, column_flags::unsigned_);
    output.emplace_back("unsigned_lt_min", "-18446744073709551616", t, column_flags::unsigned_);
    output.emplace_back("unsigned_gt_max", "18446744073709551616", t, column_flags::unsigned_);
}

void add_bit_samples(
    std::vector<error_sample>& output
)
{
    output.emplace_back("bit_string_view_too_short", "", protocol_field_type::bit, column_flags::unsigned_);
    output.emplace_back("bit_string_view_too_long", "123456789", protocol_field_type::bit, column_flags::unsigned_);
}

void add_float_samples(
    protocol_field_type t,
    string_view lt_min,
    string_view gt_max,
    std::vector<error_sample>& output
)
{
    output.emplace_back("blank", "", t);
    output.emplace_back("non_number", "abtrf", t);
    output.emplace_back("lt_min", lt_min, t);
    output.emplace_back("gt_max", gt_max, t);
    output.emplace_back("inf", "inf", t); // inf values not allowed by SQL std
    output.emplace_back("minus_inf", "-inf", t);
    output.emplace_back("nan", "nan", t); // nan values not allowed by SQL std
    output.emplace_back("minus_nan", "-nan", t);
}

void add_date_samples(std::vector<error_sample>& output)
{
    output.emplace_back("empty",            "", protocol_field_type::date);
    output.emplace_back("too_short",        "2020-05-2", protocol_field_type::date);
    output.emplace_back("too_long",         "02020-05-02", protocol_field_type::date);
    output.emplace_back("bad_delimiter",    "2020:05:02", protocol_field_type::date);
    output.emplace_back("too_many_groups",  "20-20-05-2", protocol_field_type::date);
    output.emplace_back("too_few_groups",   "2020-00005", protocol_field_type::date);
    output.emplace_back("incomplete_year",  "999-05-005", protocol_field_type::date);
    output.emplace_back("hex",              "ffff-ff-ff", protocol_field_type::date);
    output.emplace_back("null_value",       makesv("2020-05-\02"), protocol_field_type::date);
    output.emplace_back("long_year",     "10000-05-2", protocol_field_type::date);
    output.emplace_back("long_month",       "2010-005-2", protocol_field_type::date);
    output.emplace_back("long_day",         "2010-5-002", protocol_field_type::date);
    output.emplace_back("negative_year",    "-001-05-02", protocol_field_type::date);
    output.emplace_back("invalid_month",    "2010-13-02", protocol_field_type::date);
    output.emplace_back("invalid_month_max","2010-99-02", protocol_field_type::date);
    output.emplace_back("negative_month",   "2010--5-02", protocol_field_type::date);
    output.emplace_back("invalid_day",      "2010-05-32", protocol_field_type::date);
    output.emplace_back("invalid_day_max",  "2010-05-99", protocol_field_type::date);
    output.emplace_back("negative_day",     "2010-05--2", protocol_field_type::date);
}

void add_datetime_samples(
    protocol_field_type t,
    std::vector<error_sample>& output
)
{
    output.emplace_back("empty",            "", t);
    output.emplace_back("too_short_0",      "2020-05-02 23:01:0", t, 0, 0);
    output.emplace_back("too_short_1",      "2020-05-02 23:01:0.1", t, 0, 1);
    output.emplace_back("too_short_2",      "2020-05-02 23:01:00.1", t, 0, 2);
    output.emplace_back("too_short_3",      "2020-05-02 23:01:00.11", t, 0, 3);
    output.emplace_back("too_short_4",      "2020-05-02 23:01:00.111", t, 0, 4);
    output.emplace_back("too_short_5",      "2020-05-02 23:01:00.1111", t, 0, 5);
    output.emplace_back("too_short_6",      "2020-05-02 23:01:00.11111", t, 0, 6);
    output.emplace_back("too_long_0",       "2020-05-02 23:01:00.8", t, 0, 0);
    output.emplace_back("too_long_1",       "2020-05-02 23:01:00.98", t, 0, 1);
    output.emplace_back("too_long_2",       "2020-05-02 23:01:00.998", t, 0, 2);
    output.emplace_back("too_long_3",       "2020-05-02 23:01:00.9998", t, 0, 3);
    output.emplace_back("too_long_4",       "2020-05-02 23:01:00.99998", t, 0, 4);
    output.emplace_back("too_long_5",       "2020-05-02 23:01:00.999998", t, 0, 5);
    output.emplace_back("too_long_6",       "2020-05-02 23:01:00.9999998", t, 0, 6);
    output.emplace_back("no_decimals_1",    "2020-05-02 23:01:00  ", t, 0, 1);
    output.emplace_back("no_decimals_2",    "2020-05-02 23:01:00   ", t, 0, 2);
    output.emplace_back("no_decimals_3",    "2020-05-02 23:01:00     ", t, 0, 3);
    output.emplace_back("no_decimals_4",    "2020-05-02 23:01:00      ", t, 0, 4);
    output.emplace_back("no_decimals_5",    "2020-05-02 23:01:00       ", t, 0, 5);
    output.emplace_back("no_decimals_6",    "2020-05-02 23:01:00        ", t, 0, 6);
    output.emplace_back("trailing_0",       "2020-05-02 23:01:0p", t, 0, 0);
    output.emplace_back("trailing_1",       "2020-05-02 23:01:00.p", t, 0, 1);
    output.emplace_back("trailing_2",       "2020-05-02 23:01:00.1p", t, 0, 2);
    output.emplace_back("trailing_3",       "2020-05-02 23:01:00.12p", t, 0, 3);
    output.emplace_back("trailing_4",       "2020-05-02 23:01:00.123p", t, 0, 4);
    output.emplace_back("trailing_5",       "2020-05-02 23:01:00.1234p", t, 0, 5);
    output.emplace_back("trailing_6",       "2020-05-02 23:01:00.12345p", t, 0, 6);
    output.emplace_back("bad_delimiter",    "2020-05-02 23-01-00", t);
    output.emplace_back("missing_1gp_0",    "2020-05-02 23:01:  ", t);
    output.emplace_back("missing_2gp_0",    "2020-05-02 23:     ", t);
    output.emplace_back("missing_3gp_0",    "2020-05-02         ", t);
    output.emplace_back("missing_1gp_1",    "2020-05-02 23:01:.9  ", t);
    output.emplace_back("missing_2gp_1",    "2020-05-02 23:.9     ", t);
    output.emplace_back("missing_3gp_1",    "2020-05-02.9         ", t);
    output.emplace_back("invalid_year",     "10000-05-02 24:20:20.1", t, 0, 2);
    output.emplace_back("negative_year",    "-100-05-02 24:20:20", t);
    output.emplace_back("invalid_month",    "2020-13-02 24:20:20", t);
    output.emplace_back("negative_month",   "2020--5-02 24:20:20", t);
    output.emplace_back("invalid_day",      "2020-05-32 24:20:20", t);
    output.emplace_back("negative_day",     "2020-05--2 24:20:20", t);
    output.emplace_back("invalid_hour",     "2020-05-02 24:20:20", t);
    output.emplace_back("negative_hour",    "2020-05-02 -2:20:20", t);
    output.emplace_back("invalid_min",      "2020-05-02 22:60:20", t);
    output.emplace_back("negative_min",     "2020-05-02 22:-1:20", t);
    output.emplace_back("invalid_sec",      "2020-05-02 22:06:60", t);
    output.emplace_back("negative_sec",     "2020-05-02 22:06:-1", t);
    output.emplace_back("negative_micro_2", "2020-05-02 22:06:01.-1", t, 0, 2);
    output.emplace_back("negative_micro_3", "2020-05-02 22:06:01.-12", t, 0, 3);
    output.emplace_back("negative_micro_4", "2020-05-02 22:06:01.-123", t, 0, 4);
    output.emplace_back("negative_micro_5", "2020-05-02 22:06:01.-1234", t, 0, 5);
    output.emplace_back("negative_micro_6", "2020-05-02 22:06:01.-12345", t, 0, 6);
}

void add_time_samples(std::vector<error_sample>& output)
{
    output.emplace_back("empty",           "", protocol_field_type::time);
    output.emplace_back("not_numbers",     "abjkjdb67", protocol_field_type::time);
    output.emplace_back("too_short_0",     "1:20:20", protocol_field_type::time);
    output.emplace_back("too_short_1",     "1:20:20.1", protocol_field_type::time, 0, 1);
    output.emplace_back("too_short_2",     "01:20:20.1", protocol_field_type::time, 0, 2);
    output.emplace_back("too_short_3",     "01:20:20.12", protocol_field_type::time, 0, 3);
    output.emplace_back("too_short_4",     "01:20:20.123", protocol_field_type::time, 0, 4);
    output.emplace_back("too_short_5",     "01:20:20.1234", protocol_field_type::time, 0, 5);
    output.emplace_back("too_short_6",     "01:20:20.12345", protocol_field_type::time, 0, 6);
    output.emplace_back("too_long_0",      "-9999:40:40", protocol_field_type::time, 0, 0);
    output.emplace_back("too_long_1",      "-9999:40:40.1", protocol_field_type::time, 0, 1);
    output.emplace_back("too_long_2",      "-9999:40:40.12", protocol_field_type::time, 0, 2);
    output.emplace_back("too_long_3",      "-9999:40:40.123", protocol_field_type::time, 0, 3);
    output.emplace_back("too_long_4",      "-9999:40:40.1234", protocol_field_type::time, 0, 4);
    output.emplace_back("too_long_5",      "-9999:40:40.12345", protocol_field_type::time, 0, 5);
    output.emplace_back("too_long_6",      "-9999:40:40.123456", protocol_field_type::time, 0, 6);
    output.emplace_back("extra_long",      "-99999999:40:40.12345678", protocol_field_type::time, 0, 6);
    output.emplace_back("extra_long2",     "99999999999:40:40", protocol_field_type::time, 0, 6);
    output.emplace_back("decimals_0",      "01:20:20.1", protocol_field_type::time, 0, 0);
    output.emplace_back("no_decimals_1",   "01:20:20  ", protocol_field_type::time, 0, 1);
    output.emplace_back("no_decimals_2",   "01:20:20   ", protocol_field_type::time, 0, 2);
    output.emplace_back("no_decimals_3",   "01:20:20    ", protocol_field_type::time, 0, 3);
    output.emplace_back("no_decimals_4",   "01:20:20     ", protocol_field_type::time, 0, 4);
    output.emplace_back("no_decimals_5",   "01:20:20      ", protocol_field_type::time, 0, 5);
    output.emplace_back("no_decimals_6",   "01:20:20       ", protocol_field_type::time, 0, 6);
    output.emplace_back("bad_delimiter",   "01-20-20", protocol_field_type::time);
    output.emplace_back("missing_1gp_0",   "23:01:  ", protocol_field_type::time);
    output.emplace_back("missing_2gp_0",   "23:     ", protocol_field_type::time);
    output.emplace_back("missing_1gp_1",   "23:01:.9  ", protocol_field_type::time, 0, 1);
    output.emplace_back("missing_2gp_1",   "23:.9     ", protocol_field_type::time, 0, 1);
    output.emplace_back("invalid_min",     "22:60:20", protocol_field_type::time);
    output.emplace_back("negative_min",    "22:-1:20", protocol_field_type::time);
    output.emplace_back("invalid_sec",     "22:06:60", protocol_field_type::time);
    output.emplace_back("negative_sec",    "22:06:-1", protocol_field_type::time);
    output.emplace_back("invalid_micro_1", "22:06:01.99", protocol_field_type::time, 0, 1);
    output.emplace_back("invalid_micro_2", "22:06:01.999", protocol_field_type::time, 0, 2);
    output.emplace_back("invalid_micro_3", "22:06:01.9999", protocol_field_type::time, 0, 3);
    output.emplace_back("invalid_micro_4", "22:06:01.99999", protocol_field_type::time, 0, 4);
    output.emplace_back("invalid_micro_5", "22:06:01.999999", protocol_field_type::time, 0, 5);
    output.emplace_back("invalid_micro_6", "22:06:01.9999999", protocol_field_type::time, 0, 6);
    output.emplace_back("negative_micro",  "22:06:01.-1", protocol_field_type::time, 0, 2);
    output.emplace_back("lt_min",          "-900:00:00.00", protocol_field_type::time, 0, 2);
    output.emplace_back("gt_max",          "900:00:00.00", protocol_field_type::time, 0, 2);
    output.emplace_back("invalid_sign",    "x670:00:00.00", protocol_field_type::time, 0, 2);
    output.emplace_back("null_char",       makesv("20:00:\00.00"), protocol_field_type::time, 0, 2);
    output.emplace_back("trailing_0",      "22:06:01k", protocol_field_type::time, 0, 0);
    output.emplace_back("trailing_1",      "22:06:01.1k", protocol_field_type::time, 0, 1);
    output.emplace_back("trailing_2",      "22:06:01.12k", protocol_field_type::time, 0, 2);
    output.emplace_back("trailing_3",      "22:06:01.123k", protocol_field_type::time, 0, 3);
    output.emplace_back("trailing_4",      "22:06:01.1234k", protocol_field_type::time, 0, 4);
    output.emplace_back("trailing_5",      "22:06:01.12345k", protocol_field_type::time, 0, 5);
    output.emplace_back("trailing_6",      "22:06:01.123456k", protocol_field_type::time, 0, 6);
    output.emplace_back("double_sign",     "--22:06:01.123456", protocol_field_type::time, 0, 6);
}

std::vector<error_sample> make_all_samples()
{
    std::vector<error_sample> res;
    add_int_samples(protocol_field_type::tiny, res);
    add_int_samples(protocol_field_type::short_, res);
    add_int_samples(protocol_field_type::int24, res);
    add_int_samples(protocol_field_type::long_, res);
    add_int_samples(protocol_field_type::longlong, res);
    add_int_samples(protocol_field_type::year, res);
    add_bit_samples(res);
    add_float_samples(protocol_field_type::float_, "-2e90", "2e90", res);
    add_float_samples(protocol_field_type::double_, "-2e9999", "2e9999", res);
    add_date_samples(res);
    add_datetime_samples(protocol_field_type::datetime, res);
    add_datetime_samples(protocol_field_type::timestamp, res);
    add_time_samples(res);
    return res;
}

BOOST_DATA_TEST_CASE(error, data::make(make_all_samples()))
{
    auto meta = create_meta(sample.type, sample.flags, static_cast<std::uint8_t>(sample.decimals));
    auto buffer_first = reinterpret_cast<const std::uint8_t*>(sample.from.data());
    field_view actual_value;
    auto err = deserialize_text_field(sample.from, meta, buffer_first, actual_value);
    BOOST_TEST(err == sample.expected_err);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

// clang-format on

}  // namespace
