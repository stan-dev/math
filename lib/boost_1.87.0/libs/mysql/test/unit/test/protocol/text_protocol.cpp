//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/blob_view.hpp>
#include <boost/mysql/column_type.hpp>
#include <boost/mysql/date.hpp>
#include <boost/mysql/datetime.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/metadata.hpp>

#include <boost/mysql/impl/internal/protocol/impl/text_protocol.hpp>

#include <boost/test/data/monomorphic/collection.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "operators.hpp"
#include "test_common/create_basic.hpp"
#include "test_unit/create_meta.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;
using namespace boost::unit_test;

namespace {

// clang-format off

BOOST_AUTO_TEST_SUITE(test_text_protocol)

// Success cases
BOOST_AUTO_TEST_SUITE(deserialize_success)

struct success_sample
{
    std::string name;
    std::string from;
    field_view expected;
    metadata meta;

    template <class T>
    success_sample(
        std::string name,
        std::string from,
        T&& expected_value,
        metadata meta
    ) :
        name(std::move(name)),
        from(std::move(from)),
        expected(std::forward<T>(expected_value)),
        meta(std::move(meta))
    {
    }
};

std::ostream& operator<<(std::ostream& os, const success_sample& input)
{
    return os << "(input=" << input.from
              << ", type=" << input.meta.type()
              << ", name=" << input.name
              << ")";
}

void add_string_samples(std::vector<success_sample>& output)
{
    output.emplace_back("varchar_non_empty", "string", "string", create_meta(column_type::varchar));
    output.emplace_back("varchar_empty", "", "", create_meta(column_type::varchar));
    output.emplace_back("char", "", "", create_meta(column_type::char_));
    output.emplace_back("text", "value", "value", create_meta(column_type::text));
    output.emplace_back("enum", "value", "value", create_meta(column_type::enum_));
    output.emplace_back("set", "value1,value2", "value1,value2", create_meta(column_type::set));
    output.emplace_back("decimal", "1", "1", create_meta(column_type::decimal));
    output.emplace_back("json", "{}", "{}", create_meta(column_type::json));
}

void add_blob_samples(std::vector<success_sample>& output)
{
    static constexpr std::uint8_t buff [] = { 0x00, 0x01, 0x02, 0x03 };
    std::string from { 0x00, 0x01, 0x02, 0x03 };

    output.emplace_back("varbinary_non_empty", from, blob_view(buff), create_meta(column_type::varbinary));
    output.emplace_back("varbinary_empty", "", blob_view(), create_meta(column_type::varbinary));
    output.emplace_back("binary", from, blob_view(buff), create_meta(column_type::binary));
    output.emplace_back("blob", from, blob_view(buff), create_meta(column_type::blob));
    output.emplace_back("geometry", from, blob_view(buff), create_meta(column_type::geometry));

    // Anything we don't know what it is, we interpret as a blob
    output.emplace_back("unknown_protocol_type", from, blob_view(buff), create_meta(column_type::unknown));
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
    column_type type,
    std::vector<success_sample>& output
)
{
    output.emplace_back("signed", "20", std::int64_t(20), create_meta(type));
    output.emplace_back("signed_max", std::move(signed_max_s), signed_max_b, create_meta(type));
    output.emplace_back("signed_negative", "-20", std::int64_t(-20), create_meta(type));
    output.emplace_back("signed_min", std::move(signed_min_s), signed_min_b, create_meta(type));
    output.emplace_back("unsigned", "20", std::uint64_t(20), meta_builder().type(type).unsigned_flag(true).build());
    output.emplace_back("unsigned_min", "0", std::uint64_t(0), meta_builder().type(type).unsigned_flag(true).build());
    output.emplace_back("unsigned_max", std::move(unsigned_max_s),
        unsigned_max_b, meta_builder().type(type).unsigned_flag(true).build());
    output.emplace_back("unsigned_zerofill", std::move(zerofill_s),
        zerofill_b, meta_builder().type(type).unsigned_flag(true).zerofill(true).build());
}

void add_int_samples(std::vector<success_sample>& output)
{
    add_int_samples_helper(
        "127", 127,
        "-128", -128,
        "255", 255,
        "010", 10,
        column_type::tinyint,
        output
    );
    add_int_samples_helper(
        "32767", 32767,
        "-32768", -32768,
        "65535", 65535,
        "00535", 535,
        column_type::smallint,
        output
    );
    add_int_samples_helper(
        "8388607", 8388607,
        "-8388608",-8388608,
        "16777215", 16777215,
        "00007215", 7215,
        column_type::mediumint,
        output
    );
    add_int_samples_helper(
        "2147483647", 2147483647,
        "-2147483648", -2147483647 - 1, // minus is not part of literal, avoids warning
        "4294967295", 4294967295,
        "0000067295", 67295,
        column_type::int_,
        output
    );
    add_int_samples_helper(
        "9223372036854775807", 9223372036854775807,
        "-9223372036854775808", -9223372036854775807 - 1, // minus is not part of literal, avoids warning
        "18446744073709551615", 18446744073709551615ULL, // suffix needed to avoid warning
        "0000067295", 67295,
        column_type::bigint,
        output
    );

    // YEAR
    auto year_meta = meta_builder().type(column_type::year).unsigned_flag(true).build();
    output.emplace_back("regular_value", "1999", std::uint64_t(1999), year_meta);
    output.emplace_back("min", "1901", std::uint64_t(1901), year_meta);
    output.emplace_back("max", "2155", std::uint64_t(2155), year_meta);
    output.emplace_back("zero", "0000", std::uint64_t(0), year_meta);
}

// bit
void add_bit_types(std::vector<success_sample>& output)
{
    auto bit_meta = meta_builder().type(column_type::bit).unsigned_flag(true).build();
    output.emplace_back("bit_8", "\x12", std::uint64_t(0x12), bit_meta);
    output.emplace_back("bit_16", "\x12\x34", std::uint64_t(0x1234), bit_meta);
    output.emplace_back("bit_24", "\x12\x34\x56", std::uint64_t(0x123456), bit_meta);
    output.emplace_back("bit_32", "\x12\x34\x56\x78", std::uint64_t(0x12345678), bit_meta);
    output.emplace_back("bit_40", "\x12\x34\x56\x78\x9a", std::uint64_t(0x123456789a), bit_meta);
    output.emplace_back("bit_48", "\x12\x34\x56\x78\x9a\xbc", std::uint64_t(0x123456789abc), bit_meta);
    output.emplace_back("bit_56", "\x12\x34\x56\x78\x9a\xbc\xde", std::uint64_t(0x123456789abcde), bit_meta);
    output.emplace_back("bit_64", "\x12\x34\x56\x78\x9a\xbc\xde\xf0", std::uint64_t(0x123456789abcdef0), bit_meta);
}

template <class T>
void add_float_samples(
    column_type type,
    std::vector<success_sample>& output
)
{
    output.emplace_back("zero", "0", T(0.0), create_meta(type));
    output.emplace_back("integer_positive", "4", T(4.0), create_meta(type));
    output.emplace_back("integer_negative", "-5", T(-5.0), create_meta(type));
    output.emplace_back("fractional_positive", "3.147", T(3.147), create_meta(type));
    output.emplace_back("fractional_negative", "-3.147", T(-3.147), create_meta(type));
    output.emplace_back("positive_exponent_positive_integer", "3e20", T(3e20), create_meta(type));
    output.emplace_back("positive_exponent_negative_integer", "-3e20", T(-3e20), create_meta(type));
    output.emplace_back("positive_exponent_positive_fractional", "3.14e20", T(3.14e20), create_meta(type));
    output.emplace_back("positive_exponent_negative_fractional", "-3.45e20", T(-3.45e20), create_meta(type));
    output.emplace_back("negative_exponent_positive_integer", "3e-20", T(3e-20), create_meta(type));
    output.emplace_back("negative_exponent_negative_integer", "-3e-20", T(-3e-20), create_meta(type));
    output.emplace_back("negative_exponent_positive_fractional", "3.14e-20", T(3.14e-20), create_meta(type));
    output.emplace_back("negative_exponent_negative_fractional", "-3.45e-20", T(-3.45e-20), create_meta(type));
}

void add_date_samples(std::vector<success_sample>& output)
{
    output.emplace_back("regular_date", "2019-02-28", date(2019u, 2u, 28u), create_meta(column_type::date));
    output.emplace_back("leap_year", "1788-02-29", date(1788u, 2u, 29u), create_meta(column_type::date));
    output.emplace_back("min", "0000-01-01", date(0u, 1u, 1u), create_meta(column_type::date));
    output.emplace_back("max", "9999-12-31", date(9999u, 12u, 31u), create_meta(column_type::date));
    output.emplace_back("zero", "0000-00-00", date(), create_meta(column_type::date));
    output.emplace_back("zero_month", "0000-00-01", date(0u, 0u, 1u), create_meta(column_type::date));
    output.emplace_back("zero_day", "0000-01-00", date(0u, 1u, 0u), create_meta(column_type::date));
    output.emplace_back("zero_month_day_nonzero_year", "2010-00-00", date(2010u, 0u, 0u), create_meta(column_type::date));
    output.emplace_back("invalid_date", "2010-11-31", date(2010u, 11u, 31u), create_meta(column_type::date));
}

void add_datetime_samples(
    column_type type,
    std::vector<success_sample>& output
)
{
    auto meta_0decimals = meta_builder().type(type).decimals(0).build();
    output.emplace_back("0_decimals_date", "2010-02-15 00:00:00", datetime(2010u, 2u, 15u), meta_0decimals);
    output.emplace_back("0_decimals_h", "2010-02-15 02:00:00", datetime(2010u, 2u, 15u, 2u), meta_0decimals);
    output.emplace_back("0_decimals_hm", "2010-02-15 02:05:00", datetime(2010u, 2u, 15u, 2u, 5u), meta_0decimals);
    output.emplace_back("0_decimals_hms", "2010-02-15 02:05:30", datetime(2010u, 2u, 15u, 2u, 5u, 30u), meta_0decimals);
    output.emplace_back("0_decimals_min", "0000-01-01 00:00:00", datetime(0u, 1u, 1u), meta_0decimals);
    output.emplace_back("0_decimals_max", "9999-12-31 23:59:59", datetime(9999u, 12u, 31u, 23u, 59u, 59u), meta_0decimals);

    auto meta_1decimal = meta_builder().type(type).decimals(1).build();
    output.emplace_back("1_decimals_date", "2010-02-15 00:00:00.0", datetime(2010u, 2u, 15u), meta_1decimal);
    output.emplace_back("1_decimals_h", "2010-02-15 02:00:00.0", datetime(2010u, 2u, 15u, 2u), meta_1decimal);
    output.emplace_back("1_decimals_hm", "2010-02-15 02:05:00.0", datetime(2010u, 2u, 15u, 2u, 5u), meta_1decimal);
    output.emplace_back("1_decimals_hms", "2010-02-15 02:05:30.0", datetime(2010u, 2u, 15u, 2u, 5u, 30u), meta_1decimal);
    output.emplace_back("1_decimals_hmsu", "2010-02-15 02:05:30.5", datetime(2010u, 2u, 15u, 2u, 5u, 30u, 500000u), meta_1decimal);
    output.emplace_back("1_decimals_min", "0000-01-01 00:00:00.0", datetime(0u, 1u, 1u), meta_1decimal);
    output.emplace_back("1_decimals_max", "9999-12-31 23:59:59.9", datetime(9999u, 12u, 31u, 23u, 59u, 59u, 900000u), meta_1decimal);

    auto meta_2decimals = meta_builder().type(type).decimals(2).build();
    output.emplace_back("2_decimals_hms", "2010-02-15 02:05:30.00", datetime(2010u, 2u, 15u, 2u, 5u, 30u), meta_2decimals);
    output.emplace_back("2_decimals_hmsu", "2010-02-15 02:05:30.05", datetime(2010u, 2u, 15u, 2u, 5u, 30u, 50000u), meta_2decimals);
    output.emplace_back("2_decimals_min", "0000-01-01 00:00:00.00", datetime(0u, 1u, 1u), meta_2decimals);
    output.emplace_back("2_decimals_max", "9999-12-31 23:59:59.99", datetime(9999u, 12u, 31u, 23u, 59u, 59u, 990000u), meta_2decimals);

    auto meta_3decimals = meta_builder().type(type).decimals(3).build();
    output.emplace_back("3_decimals_hms", "2010-02-15 02:05:30.000", datetime(2010u, 2u, 15u, 2u, 5u, 30u), meta_3decimals);
    output.emplace_back("3_decimals_hmsu", "2010-02-15 02:05:30.420", datetime(2010u, 2u, 15u, 2u, 5u, 30u, 420000u), meta_3decimals);
    output.emplace_back("3_decimals_min", "0000-01-01 00:00:00.000", datetime(0u, 1u, 1u), meta_3decimals);
    output.emplace_back("3_decimals_max", "9999-12-31 23:59:59.999", datetime(9999u, 12u, 31u, 23u, 59u, 59u, 999000u), meta_3decimals);

    auto meta_4decimals = meta_builder().type(type).decimals(4).build();
    output.emplace_back("4_decimals_hms", "2010-02-15 02:05:30.0000", datetime(2010u, 2u, 15u, 2u, 5u, 30u), meta_4decimals);
    output.emplace_back("4_decimals_hmsu", "2010-02-15 02:05:30.4267", datetime(2010u, 2u, 15u, 2u, 5u, 30u, 426700u), meta_4decimals);
    output.emplace_back("4_decimals_min", "0000-01-01 00:00:00.0000", datetime(0u, 1u, 1u), meta_4decimals);
    output.emplace_back("4_decimals_max", "9999-12-31 23:59:59.9999", datetime(9999u, 12u, 31u, 23u, 59u, 59u, 999900u), meta_4decimals);

    auto meta_5decimals = meta_builder().type(type).decimals(5).build();
    output.emplace_back("5_decimals_hms", "2010-02-15 02:05:30.00000", datetime(2010u, 2u, 15u, 2u, 5u, 30u), meta_5decimals);
    output.emplace_back("5_decimals_hmsu", "2010-02-15 02:05:30.00239", datetime(2010u, 2u, 15u, 2u, 5u, 30u, 2390u), meta_5decimals);
    output.emplace_back("5_decimals_min", "0000-01-01 00:00:00.00000", datetime(0u, 1u, 1u), meta_5decimals);
    output.emplace_back("5_decimals_max", "9999-12-31 23:59:59.99999", datetime(9999u, 12u, 31u, 23u, 59u, 59u, 999990u), meta_5decimals);

    auto meta_6decimals = meta_builder().type(type).decimals(6).build();
    output.emplace_back("6_decimals_hms", "2010-02-15 02:05:30.000000", datetime(2010u, 2u, 15u, 2u, 5u, 30u), meta_6decimals);
    output.emplace_back("6_decimals_hmsu", "2010-02-15 02:05:30.002395", datetime(2010u, 2u, 15u, 2u, 5u, 30u, 2395u), meta_6decimals);
    output.emplace_back("6_decimals_min", "0000-01-01 00:00:00.000000", datetime(0u, 1u, 1u), meta_6decimals);
    output.emplace_back("6_decimals_max", "9999-12-31 23:59:59.999999", datetime(9999u, 12u, 31u, 23u, 59u, 59u, 999999u), meta_6decimals);

    // not a real case, we cap decimals to 6
    auto meta_7decimals = meta_builder().type(type).decimals(7).build();
    output.emplace_back("7_decimals", "2010-02-15 02:05:30.002395", datetime(2010, 2, 15, 2, 5, 30, 2395), meta_7decimals);

    // Invalid datetimes (because their date is invalid)
    output.emplace_back("0_decimals_zero", "0000-00-00 00:00:00", datetime(), meta_0decimals);
    output.emplace_back("0_decimals_invalid_date", "2010-11-31 01:10:59", datetime(2010u, 11u, 31u, 1u, 10u, 59u), meta_0decimals);
    output.emplace_back("0_decimals_zero_month", "2010-00-31 01:10:59", datetime(2010u, 0u, 31u, 1u, 10u, 59u), meta_0decimals);
    output.emplace_back("0_decimals_zero_day", "2010-11-00 01:10:59", datetime(2010u, 11u, 0u, 1u, 10u, 59u), meta_0decimals);
    output.emplace_back("0_decimals_zero_month_day", "2010-00-00 01:10:59", datetime(2010u, 0u, 0u, 1u, 10u, 59u), meta_0decimals);

    output.emplace_back("1_decimals_zero", "0000-00-00 00:00:00.0", datetime(), meta_1decimal);
    output.emplace_back("1_decimals_invalid_date", "2010-11-31 01:10:59.9", datetime(2010u, 11u, 31u, 1u, 10u, 59u, 900000u), meta_1decimal);
    output.emplace_back("1_decimals_zero_month", "2010-00-31 01:10:59.9", datetime(2010u, 0u, 31u, 1u, 10u, 59u, 900000u), meta_1decimal);
    output.emplace_back("1_decimals_zero_day", "2010-11-00 01:10:59.9", datetime(2010u, 11u, 0u, 1u, 10u, 59u, 900000u), meta_1decimal);
    output.emplace_back("1_decimals_zero_month_day", "2010-00-00 01:10:59.9", datetime(2010u, 0u, 0u, 1u, 10u, 59u, 900000u), meta_1decimal);

    output.emplace_back("2_decimals_zero", "0000-00-00 00:00:00.00", datetime(), meta_2decimals);
    output.emplace_back("2_decimals_invalid_date", "2010-11-31 01:10:59.98", datetime(2010u, 11u, 31u, 1u, 10u, 59u, 980000u), meta_2decimals);
    output.emplace_back("2_decimals_zero_month", "2010-00-31 01:10:59.98", datetime(2010u, 0u, 31u, 1u, 10u, 59u, 980000u), meta_2decimals);
    output.emplace_back("2_decimals_zero_day", "2010-11-00 01:10:59.98", datetime(2010u, 11u, 0u, 1u, 10u, 59u, 980000u), meta_2decimals);
    output.emplace_back("2_decimals_zero_month_day", "2010-00-00 01:10:59.98", datetime(2010u, 0u, 0u, 1u, 10u, 59u, 980000u), meta_2decimals);

    output.emplace_back("3_decimals_zero", "0000-00-00 00:00:00.000", datetime(), meta_3decimals);
    output.emplace_back("3_decimals_invalid_date", "2010-11-31 01:10:59.987", datetime(2010u, 11u, 31u, 1u, 10u, 59u, 987000u), meta_3decimals);
    output.emplace_back("3_decimals_zero_month", "2010-00-31 01:10:59.987", datetime(2010u, 0u, 31u, 1u, 10u, 59u, 987000u), meta_3decimals);
    output.emplace_back("3_decimals_zero_day", "2010-11-00 01:10:59.987", datetime(2010u, 11u, 0u, 1u, 10u, 59u, 987000u), meta_3decimals);
    output.emplace_back("3_decimals_zero_month_day", "2010-00-00 01:10:59.987", datetime(2010u, 0u, 0u, 1u, 10u, 59u, 987000u), meta_3decimals);

    output.emplace_back("4_decimals_zero", "0000-00-00 00:00:00.0000", datetime(), meta_4decimals);
    output.emplace_back("4_decimals_invalid_date", "2010-11-31 01:10:59.9876", datetime(2010u, 11u, 31u, 1u, 10u, 59u, 987600u), meta_4decimals);
    output.emplace_back("4_decimals_zero_month", "2010-00-31 01:10:59.9876", datetime(2010u, 0u, 31u, 1u, 10u, 59u, 987600u), meta_4decimals);
    output.emplace_back("4_decimals_zero_day", "2010-11-00 01:10:59.9876", datetime(2010u, 11u, 0u, 1u, 10u, 59u, 987600u), meta_4decimals);
    output.emplace_back("4_decimals_zero_month_day", "2010-00-00 01:10:59.9876", datetime(2010u, 0u, 0u, 1u, 10u, 59u, 987600u), meta_4decimals);

    output.emplace_back("5_decimals_zero", "0000-00-00 00:00:00.00000", datetime(), meta_5decimals);
    output.emplace_back("5_decimals_invalid_date", "2010-11-31 01:10:59.98765", datetime(2010u, 11u, 31u, 1u, 10u, 59u, 987650u), meta_5decimals);
    output.emplace_back("5_decimals_zero_month", "2010-00-31 01:10:59.98765", datetime(2010u, 0u, 31u, 1u, 10u, 59u, 987650u), meta_5decimals);
    output.emplace_back("5_decimals_zero_day", "2010-11-00 01:10:59.98765", datetime(2010u, 11u, 0u, 1u, 10u, 59u, 987650u), meta_5decimals);
    output.emplace_back("5_decimals_zero_month_day", "2010-00-00 01:10:59.98765", datetime(2010u, 0u, 0u, 1u, 10u, 59u, 987650u), meta_5decimals);

    output.emplace_back("6_decimals_zero", "0000-00-00 00:00:00.000000", datetime(), meta_6decimals);
    output.emplace_back("6_decimals_invalid_date", "2010-11-31 01:10:59.987654", datetime(2010u, 11u, 31u, 1u, 10u, 59u, 987654u), meta_6decimals);
    output.emplace_back("6_decimals_zero_month", "2010-00-31 01:10:59.987654", datetime(2010u, 0u, 31u, 1u, 10u, 59u, 987654u), meta_6decimals);
    output.emplace_back("6_decimals_zero_day", "2010-11-00 01:10:59.987654", datetime(2010u, 11u, 0u, 1u, 10u, 59u, 987654u), meta_6decimals);
    output.emplace_back("6_decimals_zero_month_day", "2010-00-00 01:10:59.987654", datetime(2010u, 0u, 0u, 1u, 10u, 59u, 987654u), meta_6decimals);
}

void add_time_samples(std::vector<success_sample>& output)
{
    auto meta_0decimals = meta_builder().type(column_type::time).decimals(0).build();
    output.emplace_back("0_decimals_positive_h", "01:00:00", maket(1, 0, 0), meta_0decimals);
    output.emplace_back("0_decimals_positive_hm", "12:03:00", maket(12, 3, 0), meta_0decimals);
    output.emplace_back("0_decimals_positive_hms", "14:51:23", maket(14, 51, 23), meta_0decimals);
    output.emplace_back("0_decimals_max", "838:59:59", maket(838, 59, 59), meta_0decimals);
    output.emplace_back("0_decimals_negative_h", "-06:00:00", -maket(6, 0, 0), meta_0decimals);
    output.emplace_back("0_decimals_negative_hm", "-12:03:00", -maket(12, 3, 0), meta_0decimals);
    output.emplace_back("0_decimals_negative_hms", "-14:51:23", -maket(14, 51, 23), meta_0decimals);
    output.emplace_back("0_decimals_min", "-838:59:59", -maket(838, 59, 59), meta_0decimals);
    output.emplace_back("0_decimals_zero", "00:00:00", maket(0, 0, 0), meta_0decimals);
    output.emplace_back("0_decimals_negative_h0", "-00:51:23", -maket(0, 51, 23), meta_0decimals);

    auto meta_1decimal = meta_builder().type(column_type::time).decimals(1).build();
    output.emplace_back("1_decimals_positive_hms", "14:51:23.0", maket(14, 51, 23), meta_1decimal);
    output.emplace_back("1_decimals_positive_hmsu", "14:51:23.5", maket(14, 51, 23, 500000), meta_1decimal);
    output.emplace_back("1_decimals_max", "838:59:58.9", maket(838, 59, 58, 900000), meta_1decimal);
    output.emplace_back("1_decimals_negative_hms", "-14:51:23.0", -maket(14, 51, 23), meta_1decimal);
    output.emplace_back("1_decimals_negative_hmsu", "-14:51:23.5", -maket(14, 51, 23, 500000), meta_1decimal);
    output.emplace_back("1_decimals_min", "-838:59:58.9", -maket(838, 59, 58, 900000), meta_1decimal);
    output.emplace_back("1_decimals_zero", "00:00:00.0", maket(0, 0, 0), meta_1decimal);
    output.emplace_back("1_decimals_negative_h0", "-00:51:23.1", -maket(0, 51, 23, 100000), meta_1decimal);

    auto meta_2decimals = meta_builder().type(column_type::time).decimals(2).build();
    output.emplace_back("2_decimals_positive_hms", "14:51:23.00", maket(14, 51, 23), meta_2decimals);
    output.emplace_back("2_decimals_positive_hmsu", "14:51:23.52", maket(14, 51, 23, 520000), meta_2decimals);
    output.emplace_back("2_decimals_max", "838:59:58.99", maket(838, 59, 58, 990000), meta_2decimals);
    output.emplace_back("2_decimals_negative_hms", "-14:51:23.00", -maket(14, 51, 23), meta_2decimals);
    output.emplace_back("2_decimals_negative_hmsu", "-14:51:23.50", -maket(14, 51, 23, 500000), meta_2decimals);
    output.emplace_back("2_decimals_min", "-838:59:58.99", -maket(838, 59, 58, 990000), meta_2decimals);
    output.emplace_back("2_decimals_zero", "00:00:00.00", maket(0, 0, 0), meta_2decimals);
    output.emplace_back("2_decimals_negative_h0", "-00:51:23.12", -maket(0, 51, 23, 120000), meta_2decimals);

    auto meta_3decimals = meta_builder().type(column_type::time).decimals(3).build();
    output.emplace_back("3_decimals_positive_hms", "14:51:23.000", maket(14, 51, 23), meta_3decimals);
    output.emplace_back("3_decimals_positive_hmsu", "14:51:23.501", maket(14, 51, 23, 501000), meta_3decimals);
    output.emplace_back("3_decimals_max", "838:59:58.999", maket(838, 59, 58, 999000), meta_3decimals);
    output.emplace_back("3_decimals_negative_hms", "-14:51:23.000", -maket(14, 51, 23), meta_3decimals);
    output.emplace_back("3_decimals_negative_hmsu", "-14:51:23.003", -maket(14, 51, 23, 3000), meta_3decimals);
    output.emplace_back("3_decimals_min", "-838:59:58.999", -maket(838, 59, 58, 999000), meta_3decimals);
    output.emplace_back("3_decimals_zero", "00:00:00.000", maket(0, 0, 0), meta_3decimals);
    output.emplace_back("3_decimals_negative_h0", "-00:51:23.123", -maket(0, 51, 23, 123000), meta_3decimals);

    auto meta_4decimals = meta_builder().type(column_type::time).decimals(4).build();
    output.emplace_back("4_decimals_positive_hms", "14:51:23.0000", maket(14, 51, 23), meta_4decimals);
    output.emplace_back("4_decimals_positive_hmsu", "14:51:23.5017", maket(14, 51, 23, 501700), meta_4decimals);
    output.emplace_back("4_decimals_max", "838:59:58.9999", maket(838, 59, 58, 999900), meta_4decimals);
    output.emplace_back("4_decimals_negative_hms", "-14:51:23.0000", -maket(14, 51, 23), meta_4decimals);
    output.emplace_back("4_decimals_negative_hmsu", "-14:51:23.0038", -maket(14, 51, 23, 3800), meta_4decimals);
    output.emplace_back("4_decimals_min", "-838:59:58.9999", -maket(838, 59, 58, 999900), meta_4decimals);
    output.emplace_back("4_decimals_zero", "00:00:00.0000", maket(0, 0, 0), meta_4decimals);
    output.emplace_back("4_decimals_negative_h0", "-00:51:23.1234", -maket(0, 51, 23, 123400), meta_4decimals);

    auto meta_5decimals = meta_builder().type(column_type::time).decimals(5).build();
    output.emplace_back("5_decimals_positive_hms", "14:51:23.00000", maket(14, 51, 23), meta_5decimals);
    output.emplace_back("5_decimals_positive_hmsu", "14:51:23.50171", maket(14, 51, 23, 501710), meta_5decimals);
    output.emplace_back("5_decimals_max", "838:59:58.99999", maket(838, 59, 58, 999990), meta_5decimals);
    output.emplace_back("5_decimals_negative_hms", "-14:51:23.00000", -maket(14, 51, 23), meta_5decimals);
    output.emplace_back("5_decimals_negative_hmsu", "-14:51:23.00009", -maket(14, 51, 23, 90), meta_5decimals);
    output.emplace_back("5_decimals_min", "-838:59:58.99999", -maket(838, 59, 58, 999990), meta_5decimals);
    output.emplace_back("5_decimals_zero", "00:00:00.00000", maket(0, 0, 0), meta_5decimals);
    output.emplace_back("5_decimals_negative_h0", "-00:51:23.12345", -maket(0, 51, 23, 123450), meta_5decimals);

    auto meta_6decimals = meta_builder().type(column_type::time).decimals(6).build();
    output.emplace_back("6_decimals_positive_hms", "14:51:23.000000", maket(14, 51, 23), meta_6decimals);
    output.emplace_back("6_decimals_positive_hmsu", "14:51:23.501717", maket(14, 51, 23, 501717), meta_6decimals);
    output.emplace_back("6_decimals_max", "838:59:58.999999", maket(838, 59, 58, 999999), meta_6decimals);
    output.emplace_back("6_decimals_negative_hms", "-14:51:23.000000", -maket(14, 51, 23), meta_6decimals);
    output.emplace_back("6_decimals_negative_hmsu", "-14:51:23.900000", -maket(14, 51, 23, 900000), meta_6decimals);
    output.emplace_back("6_decimals_min", "-838:59:58.999999", -maket(838, 59, 58, 999999), meta_6decimals);
    output.emplace_back("6_decimals_zero", "00:00:00.000000", maket(0, 0, 0), meta_6decimals);
    output.emplace_back("6_decimals_negative_h0", "-00:51:23.123456", -maket(0, 51, 23, 123456), meta_6decimals);

    // This is not a real case - we cap anything above 6 decimals to 6
    auto meta_7decimals = meta_builder().type(column_type::time).decimals(7).build();
    output.emplace_back("7_decimals", "14:51:23.501717", maket(14, 51, 23, 501717), meta_7decimals);
}

std::vector<success_sample> make_all_samples()
{
    std::vector<success_sample> res;
    add_string_samples(res);
    add_blob_samples(res);
    add_int_samples(res);
    add_bit_types(res);
    add_float_samples<float>(column_type::float_, res);
    add_float_samples<double>(column_type::double_, res);
    add_date_samples(res);
    add_datetime_samples(column_type::datetime, res);
    add_datetime_samples(column_type::timestamp, res);
    add_time_samples(res);
    return res;
}

BOOST_DATA_TEST_CASE(ok, data::make(make_all_samples()))
{
    field_view actual_value;
    auto err = detail::deserialize_text_field(sample.from, sample.meta, actual_value);

    BOOST_TEST(err == detail::deserialize_errc::ok);
    BOOST_TEST(actual_value == sample.expected);
}

BOOST_AUTO_TEST_SUITE_END()

//
// Error cases
//

BOOST_AUTO_TEST_SUITE(deserialize_error)

struct error_sample
{
    std::string name;
    string_view from;
    metadata meta;
    detail::deserialize_errc expected_err;

    error_sample(
        std::string&& name,
        string_view from,
        metadata meta,
        detail::deserialize_errc expected_err = detail::deserialize_errc::protocol_value_error
    ) :
        name(std::move(name)),
        from(from),
        meta(std::move(meta)),
        expected_err(expected_err)
    {
    }
};

std::ostream& operator<<(std::ostream& os, const error_sample& input)
{
    return os << "(input=" << input.from
              << ", type=" << input.meta.type()
              << ", name=" << input.name
              << ")";
}

void add_int_samples(
    column_type t,
    std::vector<error_sample>& output
)
{
    auto meta_signed = create_meta(t);
    output.emplace_back("signed_blank", "", meta_signed);
    output.emplace_back("signed_non_number", "abtrf", meta_signed);
    output.emplace_back("signed_hex", "0x01", meta_signed);
    output.emplace_back("signed_fractional", "1.1", meta_signed);
    output.emplace_back("signed_exp", "2e10", meta_signed);
    output.emplace_back("signed_lt_min", "-9223372036854775809", meta_signed);
    output.emplace_back("signed_gt_max", "9223372036854775808", meta_signed);

    auto meta_unsigned = meta_builder().type(t).unsigned_flag(true).build();
    output.emplace_back("unsigned_blank", "", meta_unsigned);
    output.emplace_back("unsigned_non_number", "abtrf", meta_unsigned);
    output.emplace_back("unsigned_hex", "0x01", meta_unsigned);
    output.emplace_back("unsigned_fractional", "1.1", meta_unsigned);
    output.emplace_back("unsigned_exp", "2e10", meta_unsigned);
    output.emplace_back("unsigned_lt_min", "-18446744073709551616", meta_unsigned);
    output.emplace_back("unsigned_gt_max", "18446744073709551616", meta_unsigned);
}

void add_bit_samples(
    std::vector<error_sample>& output
)
{
    auto meta = meta_builder().type(column_type::bit).unsigned_flag(true).build();
    output.emplace_back("bit_string_view_too_short", "", meta);
    output.emplace_back("bit_string_view_too_long", "123456789", meta);
}

void add_float_samples(
    column_type t,
    string_view lt_min,
    string_view gt_max,
    std::vector<error_sample>& output
)
{
    auto meta = create_meta(t);
    output.emplace_back("blank", "", meta);
    output.emplace_back("non_number", "abtrf", meta);
    output.emplace_back("lt_min", lt_min, meta);
    output.emplace_back("gt_max", gt_max, meta);
    output.emplace_back("inf", "inf", meta); // inf values not allowed by SQL std
    output.emplace_back("minus_inf", "-inf", meta);
    output.emplace_back("nan", "nan", meta); // nan values not allowed by SQL std
    output.emplace_back("minus_nan", "-nan", meta);
}

void add_date_samples(std::vector<error_sample>& output)
{
    auto meta = create_meta(column_type::date);
    output.emplace_back("empty",            "", meta);
    output.emplace_back("too_short",        "2020-05-2", meta);
    output.emplace_back("too_long",         "02020-05-02", meta);
    output.emplace_back("bad_delimiter",    "2020:05:02", meta);
    output.emplace_back("too_many_groups",  "20-20-05-2", meta);
    output.emplace_back("too_few_groups",   "2020-00005", meta);
    output.emplace_back("incomplete_year",  "999-05-005", meta);
    output.emplace_back("hex",              "ffff-ff-ff", meta);
    output.emplace_back("null_value",       makesv("2020-05-\02"), meta);
    output.emplace_back("long_year",     "10000-05-2", meta);
    output.emplace_back("long_month",       "2010-005-2", meta);
    output.emplace_back("long_day",         "2010-5-002", meta);
    output.emplace_back("negative_year",    "-001-05-02", meta);
    output.emplace_back("invalid_month",    "2010-13-02", meta);
    output.emplace_back("invalid_month_max","2010-99-02", meta);
    output.emplace_back("negative_month",   "2010--5-02", meta);
    output.emplace_back("invalid_day",      "2010-05-32", meta);
    output.emplace_back("invalid_day_max",  "2010-05-99", meta);
    output.emplace_back("negative_day",     "2010-05--2", meta);
}

void add_datetime_samples(
    column_type t,
    std::vector<error_sample>& output
)
{
    auto meta_0decimals = meta_builder().type(t).decimals(0).build();
    auto meta_1decimals = meta_builder().type(t).decimals(1).build();
    auto meta_2decimals = meta_builder().type(t).decimals(2).build();
    auto meta_3decimals = meta_builder().type(t).decimals(3).build();
    auto meta_4decimals = meta_builder().type(t).decimals(4).build();
    auto meta_5decimals = meta_builder().type(t).decimals(5).build();
    auto meta_6decimals = meta_builder().type(t).decimals(6).build();

    output.emplace_back("empty",            "", meta_0decimals);
    output.emplace_back("too_short_0",      "2020-05-02 23:01:0", meta_0decimals);
    output.emplace_back("too_short_1",      "2020-05-02 23:01:0.1", meta_1decimals);
    output.emplace_back("too_short_2",      "2020-05-02 23:01:00.1", meta_2decimals);
    output.emplace_back("too_short_3",      "2020-05-02 23:01:00.11", meta_3decimals);
    output.emplace_back("too_short_4",      "2020-05-02 23:01:00.111", meta_4decimals);
    output.emplace_back("too_short_5",      "2020-05-02 23:01:00.1111", meta_5decimals);
    output.emplace_back("too_short_6",      "2020-05-02 23:01:00.11111", meta_6decimals);
    output.emplace_back("too_long_0",       "2020-05-02 23:01:00.8", meta_0decimals);
    output.emplace_back("too_long_1",       "2020-05-02 23:01:00.98", meta_1decimals);
    output.emplace_back("too_long_2",       "2020-05-02 23:01:00.998", meta_2decimals);
    output.emplace_back("too_long_3",       "2020-05-02 23:01:00.9998", meta_3decimals);
    output.emplace_back("too_long_4",       "2020-05-02 23:01:00.99998", meta_4decimals);
    output.emplace_back("too_long_5",       "2020-05-02 23:01:00.999998", meta_5decimals);
    output.emplace_back("too_long_6",       "2020-05-02 23:01:00.9999998", meta_6decimals);
    output.emplace_back("no_decimals_1",    "2020-05-02 23:01:00  ", meta_1decimals);
    output.emplace_back("no_decimals_2",    "2020-05-02 23:01:00   ", meta_2decimals);
    output.emplace_back("no_decimals_3",    "2020-05-02 23:01:00     ", meta_3decimals);
    output.emplace_back("no_decimals_4",    "2020-05-02 23:01:00      ", meta_4decimals);
    output.emplace_back("no_decimals_5",    "2020-05-02 23:01:00       ", meta_5decimals);
    output.emplace_back("no_decimals_6",    "2020-05-02 23:01:00        ", meta_6decimals);
    output.emplace_back("trailing_0",       "2020-05-02 23:01:0p", meta_0decimals);
    output.emplace_back("trailing_1",       "2020-05-02 23:01:00.p", meta_1decimals);
    output.emplace_back("trailing_2",       "2020-05-02 23:01:00.1p", meta_2decimals);
    output.emplace_back("trailing_3",       "2020-05-02 23:01:00.12p", meta_3decimals);
    output.emplace_back("trailing_4",       "2020-05-02 23:01:00.123p", meta_4decimals);
    output.emplace_back("trailing_5",       "2020-05-02 23:01:00.1234p", meta_5decimals);
    output.emplace_back("trailing_6",       "2020-05-02 23:01:00.12345p", meta_6decimals);
    output.emplace_back("bad_delimiter",    "2020-05-02 23-01-00", meta_0decimals);
    output.emplace_back("missing_1gp_0",    "2020-05-02 23:01:  ", meta_0decimals);
    output.emplace_back("missing_2gp_0",    "2020-05-02 23:     ", meta_0decimals);
    output.emplace_back("missing_3gp_0",    "2020-05-02         ", meta_0decimals);
    output.emplace_back("missing_1gp_1",    "2020-05-02 23:01:.9  ", meta_0decimals);
    output.emplace_back("missing_2gp_1",    "2020-05-02 23:.9     ", meta_0decimals);
    output.emplace_back("missing_3gp_1",    "2020-05-02.9         ", meta_0decimals);
    output.emplace_back("invalid_year",     "10000-05-02 24:20:20.1", meta_2decimals);
    output.emplace_back("negative_year",    "-100-05-02 24:20:20", meta_0decimals);
    output.emplace_back("invalid_month",    "2020-13-02 24:20:20", meta_0decimals);
    output.emplace_back("negative_month",   "2020--5-02 24:20:20", meta_0decimals);
    output.emplace_back("invalid_day",      "2020-05-32 24:20:20", meta_0decimals);
    output.emplace_back("negative_day",     "2020-05--2 24:20:20", meta_0decimals);
    output.emplace_back("invalid_hour",     "2020-05-02 24:20:20", meta_0decimals);
    output.emplace_back("negative_hour",    "2020-05-02 -2:20:20", meta_0decimals);
    output.emplace_back("invalid_min",      "2020-05-02 22:60:20", meta_0decimals);
    output.emplace_back("negative_min",     "2020-05-02 22:-1:20", meta_0decimals);
    output.emplace_back("invalid_sec",      "2020-05-02 22:06:60", meta_0decimals);
    output.emplace_back("negative_sec",     "2020-05-02 22:06:-1", meta_0decimals);
    output.emplace_back("negative_micro_2", "2020-05-02 22:06:01.-1", meta_2decimals);
    output.emplace_back("negative_micro_3", "2020-05-02 22:06:01.-12", meta_3decimals);
    output.emplace_back("negative_micro_4", "2020-05-02 22:06:01.-123", meta_4decimals);
    output.emplace_back("negative_micro_5", "2020-05-02 22:06:01.-1234", meta_5decimals);
    output.emplace_back("negative_micro_6", "2020-05-02 22:06:01.-12345", meta_6decimals);
}

void add_time_samples(std::vector<error_sample>& output)
{
    auto meta_0decimals = meta_builder().type(column_type::time).decimals(0).build();
    auto meta_1decimals = meta_builder().type(column_type::time).decimals(1).build();
    auto meta_2decimals = meta_builder().type(column_type::time).decimals(2).build();
    auto meta_3decimals = meta_builder().type(column_type::time).decimals(3).build();
    auto meta_4decimals = meta_builder().type(column_type::time).decimals(4).build();
    auto meta_5decimals = meta_builder().type(column_type::time).decimals(5).build();
    auto meta_6decimals = meta_builder().type(column_type::time).decimals(6).build();

    output.emplace_back("empty",           "", meta_0decimals);
    output.emplace_back("not_numbers",     "abjkjdb67", meta_0decimals);
    output.emplace_back("too_short_0",     "1:20:20", meta_0decimals);
    output.emplace_back("too_short_1",     "1:20:20.1", meta_1decimals);
    output.emplace_back("too_short_2",     "01:20:20.1", meta_2decimals);
    output.emplace_back("too_short_3",     "01:20:20.12", meta_3decimals);
    output.emplace_back("too_short_4",     "01:20:20.123", meta_4decimals);
    output.emplace_back("too_short_5",     "01:20:20.1234", meta_5decimals);
    output.emplace_back("too_short_6",     "01:20:20.12345", meta_6decimals);
    output.emplace_back("too_long_0",      "-9999:40:40", meta_0decimals);
    output.emplace_back("too_long_1",      "-9999:40:40.1", meta_1decimals);
    output.emplace_back("too_long_2",      "-9999:40:40.12", meta_2decimals);
    output.emplace_back("too_long_3",      "-9999:40:40.123", meta_3decimals);
    output.emplace_back("too_long_4",      "-9999:40:40.1234", meta_4decimals);
    output.emplace_back("too_long_5",      "-9999:40:40.12345", meta_5decimals);
    output.emplace_back("too_long_6",      "-9999:40:40.123456", meta_6decimals);
    output.emplace_back("extra_long",      "-99999999:40:40.12345678", meta_6decimals);
    output.emplace_back("extra_long2",     "99999999999:40:40", meta_6decimals);
    output.emplace_back("decimals_0",      "01:20:20.1", meta_0decimals);
    output.emplace_back("no_decimals_1",   "01:20:20  ", meta_1decimals);
    output.emplace_back("no_decimals_2",   "01:20:20   ", meta_2decimals);
    output.emplace_back("no_decimals_3",   "01:20:20    ", meta_3decimals);
    output.emplace_back("no_decimals_4",   "01:20:20     ", meta_4decimals);
    output.emplace_back("no_decimals_5",   "01:20:20      ", meta_5decimals);
    output.emplace_back("no_decimals_6",   "01:20:20       ", meta_6decimals);
    output.emplace_back("bad_delimiter",   "01-20-20", meta_0decimals);
    output.emplace_back("missing_1gp_0",   "23:01:  ", meta_0decimals);
    output.emplace_back("missing_2gp_0",   "23:     ", meta_0decimals);
    output.emplace_back("missing_1gp_1",   "23:01:.9  ", meta_1decimals);
    output.emplace_back("missing_2gp_1",   "23:.9     ", meta_1decimals);
    output.emplace_back("invalid_min",     "22:60:20", meta_0decimals);
    output.emplace_back("negative_min",    "22:-1:20", meta_0decimals);
    output.emplace_back("invalid_sec",     "22:06:60", meta_0decimals);
    output.emplace_back("negative_sec",    "22:06:-1", meta_0decimals);
    output.emplace_back("invalid_micro_1", "22:06:01.99", meta_1decimals);
    output.emplace_back("invalid_micro_2", "22:06:01.999", meta_2decimals);
    output.emplace_back("invalid_micro_3", "22:06:01.9999", meta_3decimals);
    output.emplace_back("invalid_micro_4", "22:06:01.99999", meta_4decimals);
    output.emplace_back("invalid_micro_5", "22:06:01.999999", meta_5decimals);
    output.emplace_back("invalid_micro_6", "22:06:01.9999999", meta_6decimals);
    output.emplace_back("negative_micro",  "22:06:01.-1", meta_2decimals);
    output.emplace_back("lt_min",          "-900:00:00.00", meta_2decimals);
    output.emplace_back("gt_max",          "900:00:00.00", meta_2decimals);
    output.emplace_back("invalid_sign",    "x670:00:00.00", meta_2decimals);
    output.emplace_back("null_char",       makesv("20:00:\00.00"), meta_2decimals);
    output.emplace_back("trailing_0",      "22:06:01k", meta_0decimals);
    output.emplace_back("trailing_1",      "22:06:01.1k", meta_1decimals);
    output.emplace_back("trailing_2",      "22:06:01.12k", meta_2decimals);
    output.emplace_back("trailing_3",      "22:06:01.123k", meta_3decimals);
    output.emplace_back("trailing_4",      "22:06:01.1234k", meta_4decimals);
    output.emplace_back("trailing_5",      "22:06:01.12345k", meta_5decimals);
    output.emplace_back("trailing_6",      "22:06:01.123456k", meta_6decimals);
    output.emplace_back("double_sign",     "--22:06:01.123456", meta_6decimals);
}

std::vector<error_sample> make_all_samples()
{
    std::vector<error_sample> res;
    add_int_samples(column_type::tinyint, res);
    add_int_samples(column_type::smallint, res);
    add_int_samples(column_type::mediumint, res);
    add_int_samples(column_type::int_, res);
    add_int_samples(column_type::bigint, res);
    add_int_samples(column_type::year, res);
    add_bit_samples(res);
    add_float_samples(column_type::float_, "-2e90", "2e90", res);
    add_float_samples(column_type::double_, "-2e9999", "2e9999", res);
    add_date_samples(res);
    add_datetime_samples(column_type::datetime, res);
    add_datetime_samples(column_type::timestamp, res);
    add_time_samples(res);
    return res;
}

BOOST_DATA_TEST_CASE(error, data::make(make_all_samples()))
{
    field_view actual_value;
    auto err = detail::deserialize_text_field(sample.from, sample.meta, actual_value);
    
    BOOST_TEST(err == sample.expected_err);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

// clang-format on

}  // namespace
