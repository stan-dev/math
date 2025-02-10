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

#include <boost/mysql/impl/internal/protocol/impl/binary_protocol.hpp>
#include <boost/mysql/impl/internal/protocol/impl/serialization_context.hpp>

#include <boost/test/unit_test.hpp>

#include <cstddef>

#include "operators.hpp"
#include "serialization_test.hpp"
#include "test_common/create_basic.hpp"
#include "test_unit/create_meta.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;
using detail::deserialize_errc;

namespace {

BOOST_AUTO_TEST_SUITE(test_binary_protocol)

// Adapt field_view to be Serializable
struct field_view_adaptor
{
    field_view fv;

    void serialize(detail::serialization_context& ctx) { detail::serialize_binary_field(ctx, fv); }
};

BOOST_AUTO_TEST_CASE(serialize)
{
    std::uint8_t blob_buffer[] = {0x70, 0x00, 0x01};
    struct
    {
        const char* name;
        field_view value;
        std::vector<std::uint8_t> serialized;
    } test_cases[]{
        // Strings and ints: extensive testing already done. Ensure
        // we call the right function
        {"string", field_view("abc"), {0x03, 0x61, 0x62, 0x63}},
        {"blob", field_view(blob_view(blob_buffer)), {0x03, 0x70, 0x00, 0x01}},
        {"uint64",
         field_view(std::uint64_t(0xf8f9fafbfcfdfeff)),
         {0xff, 0xfe, 0xfd, 0xfc, 0xfb, 0xfa, 0xf9, 0xf8}},
        {"int64",
         field_view(std::int64_t(-0x0706050403020101)),
         {0xff, 0xfe, 0xfd, 0xfc, 0xfb, 0xfa, 0xf9, 0xf8}},

        {"float_fractional_negative", field_view(-4.2f), {0x66, 0x66, 0x86, 0xc0}},
        {"float_fractional_positive", field_view(4.2f), {0x66, 0x66, 0x86, 0x40}},
        {"float_positive_exp_positive_fractional", field_view(3.14e20f), {0x01, 0x2d, 0x88, 0x61}},
        {"float_zero", field_view(0.0f), {0x00, 0x00, 0x00, 0x00}},

        {"double_fractional_negative", field_view(-4.2), {0xcd, 0xcc, 0xcc, 0xcc, 0xcc, 0xcc, 0x10, 0xc0}},
        {"double_fractional_positive", field_view(4.2), {0xcd, 0xcc, 0xcc, 0xcc, 0xcc, 0xcc, 0x10, 0x40}},
        {"double_positive_exp_positive_fractional",
         field_view(3.14e200),
         {0xce, 0x46, 0x3c, 0x76, 0x9c, 0x68, 0x90, 0x69}},
        {"double_zero", field_view(0.0), {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}},

        {"date_regular", field_view(date(2010u, 3u, 28u)), {0x04, 0xda, 0x07, 0x03, 0x1c}},
        {"date_min", field_view(date(1000u, 1u, 1u)), {0x04, 0xe8, 0x03, 0x01, 0x01}},
        {"date_mysqlmax", field_view(date(9999u, 12u, 31u)), {0x04, 0x0f, 0x27, 0x0c, 0x1f}},
        {"date_max", field_view(date(0xffff, 0xff, 0xff)), {0x04, 0xff, 0xff, 0xff, 0xff}},
        {"date_zero", field_view(date()), {0x04, 0x00, 0x00, 0x00, 0x00}},

        {"datetime_regular",
         field_view(datetime(2010u, 1u, 1u, 23u, 1u, 59u, 967510u)),
         {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00}},
        {"datetime_min",
         field_view(datetime(0, 1, 1)),
         {0x0b, 0x00, 0x00, 0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}},
        {"datetime_mysqlmax",
         field_view(datetime(9999, 12, 31, 23, 59, 59, 999999)),
         {0x0b, 0x0f, 0x27, 0x0c, 0x1f, 0x17, 0x3b, 0x3b, 0x3f, 0x42, 0x0f, 0x00}},
        {"datetime_max",
         field_view(datetime(0xffff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xffffffff)),
         {0x0b, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff}},
        {"datetime_zero",
         field_view(datetime()),
         {0x0b, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}},

        {"time_positive_u",
         field_view(maket(0, 0, 0, 321000)),
         {0x0c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xe8, 0xe5, 0x04, 0x00}},
        {"time_positive_hmsu",
         field_view(maket(838, 59, 58, 999000)),
         {0x0c, 0x00, 0x22, 0x00, 0x00, 0x00, 0x16, 0x3b, 0x3a, 0x58, 0x3e, 0x0f, 0x00}},
        {"time_negative_u",
         field_view(-maket(0, 0, 0, 321000)),
         {0x0c, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xe8, 0xe5, 0x04, 0x00}},
        {"time_negative_hmsu",
         field_view(-maket(838, 59, 58, 999000)),
         {0x0c, 0x01, 0x22, 0x00, 0x00, 0x00, 0x16, 0x3b, 0x3a, 0x58, 0x3e, 0x0f, 0x00}},

        // NULL is transmitted as the NULL bitmap, so nothing is expected as output
        {"null", field_view(), {}},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name) { do_serialize_test(field_view_adaptor{tc.value}, tc.serialized); }
    }
}

BOOST_AUTO_TEST_SUITE(deserialize_success)

struct success_sample
{
    std::string name;
    deserialization_buffer from;
    field_view expected;
    metadata meta;

    template <class T>
    success_sample(std::string name, std::vector<std::uint8_t> from, T&& expected_value, metadata meta)
        : name(std::move(name)),
          from(std::move(from)),
          expected(std::forward<T>(expected_value)),
          meta(std::move(meta))
    {
    }
};

void add_string_samples(std::vector<success_sample>& output)
{
    output.push_back(
        success_sample("varchar", {0x04, 0x74, 0x65, 0x73, 0x74}, "test", create_meta(column_type::varchar))
    );
    output.push_back(
        success_sample("char", {0x04, 0x74, 0x65, 0x73, 0x74}, "test", create_meta(column_type::char_))
    );
    output.push_back(
        success_sample("text", {0x04, 0x74, 0x65, 0x73, 0x74}, "test", create_meta(column_type::text))
    );
    output.push_back(
        success_sample("enum", {0x04, 0x74, 0x65, 0x73, 0x74}, "test", create_meta(column_type::enum_))
    );
    output.push_back(
        success_sample("set", {0x04, 0x74, 0x65, 0x73, 0x74}, "test", create_meta(column_type::set))
    );
    output.push_back(success_sample("decimal", {0x02, 0x31, 0x30}, "10", create_meta(column_type::decimal)));
    output.push_back(success_sample("json", {0x02, 0x7b, 0x7d}, "{}", create_meta(column_type::json)));
}

void add_blob_samples(std::vector<success_sample>& output)
{
    static constexpr std::uint8_t buff[] = {0x01, 0x00, 0x73, 0x74};

    output.push_back(success_sample(
        "varbinary",
        {0x04, 0x01, 0x00, 0x73, 0x74},
        blob_view(buff),
        create_meta(column_type::varbinary)
    ));
    output.push_back(success_sample(
        "binary",
        {0x04, 0x01, 0x00, 0x73, 0x74},
        blob_view(buff),
        create_meta(column_type::binary)
    ));
    output.push_back(success_sample(
        "blob",
        {0x04, 0x01, 0x00, 0x73, 0x74},
        blob_view(buff),
        create_meta(column_type::blob)
    ));
    output.push_back(success_sample(
        "geometry",
        {0x04, 0x01, 0x00, 0x73, 0x74},
        blob_view(buff),
        create_meta(column_type::geometry)
    ));

    // Anything we don't know what it is, we interpret as a blob
    output.push_back(success_sample(
        "unknown_protocol_type",
        {0x04, 0x01, 0x00, 0x73, 0x74},
        blob_view(buff),
        create_meta(column_type::unknown)
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
        meta_builder().type(column_type::tinyint).unsigned_flag(true).build()
    ));
    output.push_back(
        success_sample("tinyint_signed", {0xec}, std::int64_t(-20), create_meta(column_type::tinyint))
    );

    output.push_back(success_sample(
        "smallint_unsigned",
        {0x14, 0x00},
        std::uint64_t(20),
        meta_builder().type(column_type::smallint).unsigned_flag(true).build()
    ));
    output.push_back(
        success_sample("smallint_signed", {0xec, 0xff}, std::int64_t(-20), create_meta(column_type::smallint))
    );

    output.push_back(success_sample(
        "mediumint_unsigned",
        {0x14, 0x00, 0x00, 0x00},
        std::uint64_t(20),
        meta_builder().type(column_type::mediumint).unsigned_flag(true).build()
    ));
    output.push_back(success_sample(
        "mediumint_signed",
        {0xec, 0xff, 0xff, 0xff},
        std::int64_t(-20),
        create_meta(column_type::mediumint)
    ));

    output.push_back(success_sample(
        "int_unsigned",
        {0x14, 0x00, 0x00, 0x00},
        std::uint64_t(20),
        meta_builder().type(column_type::int_).unsigned_flag(true).build()
    ));
    output.push_back(success_sample(
        "int_signed",
        {0xec, 0xff, 0xff, 0xff},
        std::int64_t(-20),
        create_meta(column_type::int_)
    ));

    output.push_back(success_sample(
        "bigint_unsigned",
        {0x14, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
        std::uint64_t(20),
        meta_builder().type(column_type::bigint).unsigned_flag(true).build()
    ));
    output.push_back(success_sample(
        "bigint_signed",
        {0xec, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff},
        std::int64_t(-20),
        create_meta(column_type::bigint)
    ));

    output.push_back(success_sample(
        "year",
        {0xe3, 0x07},
        std::uint64_t(2019),
        meta_builder().type(column_type::year).unsigned_flag(true).build()
    ));
}

// bit
void add_bit_types(std::vector<success_sample>& output)
{
    auto meta = meta_builder().type(column_type::bit).unsigned_flag(true).build();

    output.push_back(success_sample("bit_8", {0x01, 0x12}, std::uint64_t(0x12), meta));
    output.push_back(success_sample("bit_16", {0x02, 0x12, 0x34}, std::uint64_t(0x1234), meta));
    output.push_back(success_sample("bit_24", {0x03, 0x12, 0x34, 0x56}, std::uint64_t(0x123456), meta));
    output.push_back(success_sample("bit_32", {0x04, 0x12, 0x34, 0x56, 0x78}, std::uint64_t(0x12345678), meta)
    );
    output.push_back(
        success_sample("bit_40", {0x05, 0x12, 0x34, 0x56, 0x78, 0x9a}, std::uint64_t(0x123456789a), meta)
    );
    output.push_back(success_sample(
        "bit_48",
        {0x06, 0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc},
        std::uint64_t(0x123456789abc),
        meta
    ));
    output.push_back(success_sample(
        "bit_56",
        {0x07, 0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde},
        std::uint64_t(0x123456789abcde),
        meta
    ));
    output.push_back(success_sample(
        "bit_64",
        {0x08, 0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xf0},
        std::uint64_t(0x123456789abcdef0),
        meta
    ));
}

void add_float_samples(std::vector<success_sample>& output)
{
    auto meta = create_meta(column_type::float_);
    output.push_back(success_sample("fractional_negative", {0x66, 0x66, 0x86, 0xc0}, -4.2f, meta));
    output.push_back(success_sample("fractional_positive", {0x66, 0x66, 0x86, 0x40}, 4.2f, meta));
    output.push_back(
        success_sample("positive_exp_positive_fractional", {0x01, 0x2d, 0x88, 0x61}, 3.14e20f, meta)
    );
    output.push_back(success_sample("zero", {0x00, 0x00, 0x00, 0x00}, 0.0f, meta));
}

void add_double_samples(std::vector<success_sample>& output)
{
    auto meta = create_meta(column_type::double_);
    output.push_back(
        success_sample("fractional_negative", {0xcd, 0xcc, 0xcc, 0xcc, 0xcc, 0xcc, 0x10, 0xc0}, -4.2, meta)
    );
    output.push_back(
        success_sample("fractional_positive", {0xcd, 0xcc, 0xcc, 0xcc, 0xcc, 0xcc, 0x10, 0x40}, 4.2, meta)
    );
    output.push_back(success_sample(
        "positive_exp_positive_fractional",
        {0xce, 0x46, 0x3c, 0x76, 0x9c, 0x68, 0x90, 0x69},
        3.14e200,
        meta
    ));
    output.push_back(success_sample("zero", {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, 0.0, meta));
}

void add_date_samples(std::vector<success_sample>& output)
{
    auto meta = create_meta(column_type::date);
    output.push_back(success_sample("regular", {0x04, 0xda, 0x07, 0x03, 0x1c}, date(2010u, 3u, 28u), meta));
    output.push_back(success_sample("min", {0x04, 0x00, 0x00, 0x01, 0x01}, date(0u, 1u, 1u), meta));
    output.push_back(success_sample("max", {0x04, 0x0f, 0x27, 0x0c, 0x1f}, date(9999u, 12u, 31u), meta));
    output.push_back(success_sample("empty", {0x00}, date(), meta));
    output.push_back(success_sample("zero", {0x04, 0x00, 0x00, 0x00, 0x00}, date(), meta));
    output.push_back(success_sample("zero_month", {0x04, 0xda, 0x07, 0x00, 0x01}, date(2010u, 0u, 1u), meta));
    output.push_back(success_sample("zero_day", {0x04, 0xda, 0x07, 0x01, 0x00}, date(2010u, 1u, 0u), meta));
    output.push_back(
        success_sample("zero_month_day", {0x04, 0xda, 0x07, 0x00, 0x00}, date(2010u, 0u, 0u), meta)
    );
    output.push_back(
        success_sample("invalid_date", {0x04, 0xda, 0x07, 0x0b, 0x1f}, date(2010u, 11u, 31u), meta)
    );
}

void add_datetime_samples(column_type type, std::vector<success_sample>& output)
{
    auto meta = create_meta(type);
    output.push_back(
        success_sample("only_date", {0x04, 0xda, 0x07, 0x01, 0x01}, datetime(2010u, 1u, 1u), meta)
    );
    output.push_back(success_sample(
        "date_h",
        {0x07, 0xda, 0x07, 0x01, 0x01, 0x14, 0x00, 0x00},
        datetime(2010u, 1u, 1u, 20u, 0u, 0u, 0u),
        meta
    ));
    output.push_back(success_sample(
        "date_m",
        {0x07, 0xda, 0x07, 0x01, 0x01, 0x00, 0x01, 0x00},
        datetime(2010u, 1u, 1u, 0u, 1u, 0u, 0u),
        meta
    ));
    output.push_back(success_sample(
        "date_hm",
        {0x07, 0xda, 0x07, 0x01, 0x01, 0x03, 0x02, 0x00},
        datetime(2010u, 1u, 1u, 3u, 2u, 0u, 0u),
        meta
    ));
    output.push_back(success_sample(
        "date_s",
        {0x07, 0xda, 0x07, 0x01, 0x01, 0x00, 0x00, 0x01},
        datetime(2010u, 1u, 1u, 0u, 0u, 1u, 0u),
        meta
    ));
    output.push_back(success_sample(
        "date_ms",
        {0x07, 0xda, 0x07, 0x01, 0x01, 0x00, 0x3b, 0x01},
        datetime(2010u, 1u, 1u, 0u, 59u, 1u, 0u),
        meta
    ));
    output.push_back(success_sample(
        "date_hs",
        {0x07, 0xda, 0x07, 0x01, 0x01, 0x05, 0x00, 0x01},
        datetime(2010u, 1u, 1u, 5u, 0u, 1u, 0u),
        meta
    ));
    output.push_back(success_sample(
        "date_hms",
        {0x07, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 0x3b},
        datetime(2010u, 1u, 1u, 23u, 1u, 59u, 0u),
        meta
    ));
    output.push_back(success_sample(
        "date_u",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x00, 0x00, 0x00, 0x78, 0xd4, 0x03, 0x00},
        datetime(2010u, 1u, 1u, 0u, 0u, 0u, 251000u),
        meta
    ));
    output.push_back(success_sample(
        "date_hu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x00, 0x00, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 1u, 1u, 23u, 0u, 0u, 967510u),
        meta
    ));
    output.push_back(success_sample(
        "date_mu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x00, 0x01, 0x00, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 1u, 1u, 0u, 1u, 0u, 967510u),
        meta
    ));
    output.push_back(success_sample(
        "date_hmu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 0x00, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 1u, 1u, 23u, 1u, 0u, 967510u),
        meta
    ));
    output.push_back(success_sample(
        "date_su",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x00, 0x00, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 1u, 1u, 0u, 0u, 59u, 967510u),
        meta
    ));
    output.push_back(success_sample(
        "date_msu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x00, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 1u, 1u, 0u, 1u, 59u, 967510u),
        meta
    ));
    output.push_back(success_sample(
        "date_hsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x00, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 1u, 1u, 23u, 0u, 59u, 967510u),
        meta
    ));
    output.push_back(success_sample(
        "date_hmsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 1u, 1u, 23u, 1u, 59u, 967510u),
        meta
    ));

    // Invalid datetimes (because their date is invalid)
    output.push_back(success_sample("empty", {0x00}, datetime(), meta));
    output.push_back(success_sample("only_date_zeros", {0x04, 0x00, 0x00, 0x00, 0x00}, datetime(), meta));
    output.push_back(success_sample(
        "only_date_invalid_date",
        {0x04, 0xda, 0x07, 0x0b, 0x1f},
        datetime(2010u, 11u, 31u),
        meta
    ));
    output.push_back(
        success_sample("only_date_zero_month", {0x04, 0xda, 0x07, 0x00, 0x01}, datetime(2010u, 0u, 1u), meta)
    );
    output.push_back(
        success_sample("only_date_zero_day", {0x04, 0xda, 0x07, 0x01, 0x00}, datetime(2010u, 1u, 0u), meta)
    );
    output.push_back(success_sample(
        "only_date_zero_month_day",
        {0x04, 0xda, 0x07, 0x00, 0x00},
        datetime(2010u, 0u, 0u),
        meta
    ));

    output.push_back(
        success_sample("date_hms_zeros", {0x07, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, datetime(), meta)
    );
    output.push_back(success_sample(
        "date_hms_invalid_date",
        {0x07, 0xda, 0x07, 0x0b, 0x1f, 0x17, 0x01, 0x3b},
        datetime(2010u, 11u, 31u, 23u, 1u, 59u),
        meta
    ));
    output.push_back(success_sample(
        "date_hms_zero_month",
        {0x07, 0xda, 0x07, 0x00, 0x01, 0x17, 0x01, 0x3b},
        datetime(2010u, 0u, 1u, 23u, 1u, 59u),
        meta
    ));
    output.push_back(success_sample(
        "date_hms_zero_day",
        {0x07, 0xda, 0x07, 0x01, 0x00, 0x17, 0x01, 0x3b},
        datetime(2010u, 1u, 0u, 23u, 1u, 59u),
        meta
    ));
    output.push_back(success_sample(
        "date_hms_zero_month_day",
        {0x07, 0xda, 0x07, 0x00, 0x00, 0x17, 0x01, 0x3b},
        datetime(2010u, 0u, 0u, 23u, 1u, 59u),
        meta
    ));

    output.push_back(success_sample(
        "date_hmsu_zeros",
        {0x0b, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
        datetime(),
        meta
    ));
    output.push_back(success_sample(
        "date_hmsu_invalid_date",
        {0x0b, 0xda, 0x07, 0x0b, 0x1f, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 11u, 31u, 23u, 1u, 59u, 967510u),
        meta
    ));
    output.push_back(success_sample(
        "date_hmsu_zero_month",
        {0x0b, 0xda, 0x07, 0x00, 0x01, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 0u, 1u, 23u, 1u, 59u, 967510u),
        meta
    ));
    output.push_back(success_sample(
        "date_hmsu_zero_day",
        {0x0b, 0xda, 0x07, 0x01, 0x00, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 1u, 0u, 23u, 1u, 59u, 967510u),
        meta
    ));
    output.push_back(success_sample(
        "date_hmsu_zero_month_day",
        {0x0b, 0xda, 0x07, 0x00, 0x00, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        datetime(2010u, 0u, 0u, 23u, 1u, 59u, 967510u),
        meta
    ));
}

void add_time_samples(std::vector<success_sample>& output)
{
    auto meta = create_meta(column_type::time);
    output.push_back(success_sample("zero", {0x00}, maket(0, 0, 0), meta));
    output.push_back(success_sample(
        "positive_d",
        {0x08, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
        maket(48, 0, 0, 0),
        meta
    ));
    output.push_back(success_sample(
        "positive_h",
        {0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x15, 0x00, 0x00},
        maket(21, 0, 0, 0),
        meta
    ));
    output.push_back(success_sample(
        "positive_m",
        {0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x28, 0x00},
        maket(0, 40, 0),
        meta
    ));
    output.push_back(success_sample(
        "positive_s",
        {0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x15},
        maket(0, 0, 21),
        meta
    ));
    output.push_back(success_sample(
        "positive_u",
        {0x0c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xe8, 0xe5, 0x04, 0x00},
        maket(0, 0, 0, 321000),
        meta
    ));
    output.push_back(success_sample(
        "positive_hmsu",
        {0x0c, 0x00, 0x22, 0x00, 0x00, 0x00, 0x16, 0x3b, 0x3a, 0x58, 0x3e, 0x0f, 0x00},
        maket(838, 59, 58, 999000),
        meta
    ));
    output.push_back(success_sample(
        "negative_d",
        {0x08, 0x01, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
        -maket(48, 0, 0, 0),
        meta
    ));
    output.push_back(success_sample(
        "negative_h",
        {0x08, 0x01, 0x00, 0x00, 0x00, 0x00, 0x15, 0x00, 0x00},
        -maket(21, 0, 0, 0),
        meta
    ));
    output.push_back(success_sample(
        "negative_m",
        {0x08, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x28, 0x00},
        -maket(0, 40, 0),
        meta
    ));
    output.push_back(success_sample(
        "negative_s",
        {0x08, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x15},
        -maket(0, 0, 21),
        meta
    ));
    output.push_back(success_sample(
        "negative_u",
        {0x0c, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xe8, 0xe5, 0x04, 0x00},
        -maket(0, 0, 0, 321000),
        meta
    ));
    output.push_back(success_sample(
        "negative_hmsu",
        {0x0c, 0x01, 0x22, 0x00, 0x00, 0x00, 0x16, 0x3b, 0x3a, 0x58, 0x3e, 0x0f, 0x00},
        -maket(838, 59, 58, 999000),
        meta
    ));
    output.push_back(success_sample(
        "negative_sign_not_one",
        {0x0c, 0x03, 0x22, 0x00, 0x00, 0x00, 0x16, 0x3b, 0x3a, 0x58, 0x3e, 0x0f, 0x00},
        -maket(838, 59, 58, 999000),
        meta
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
    add_datetime_samples(column_type::datetime, res);
    add_datetime_samples(column_type::timestamp, res);
    add_time_samples(res);
    return res;
}

BOOST_AUTO_TEST_CASE(success)
{
    for (const auto& tc : make_all_samples())
    {
        BOOST_TEST_CONTEXT("type=" << tc.meta.type() << ", name=" << tc.name)
        {
            const auto& buffer = tc.from;
            detail::deserialization_context ctx(buffer);

            field_view actual_value;
            auto err = deserialize_binary_field(ctx, tc.meta, actual_value);

            BOOST_TEST(err == deserialize_errc::ok);
            BOOST_TEST(actual_value == tc.expected);
            BOOST_TEST(ctx.first() == buffer.data() + buffer.size());  // all bytes consumed
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(deserialize_error)

struct error_sample
{
    std::string name;
    deserialization_buffer from;
    metadata meta;
    deserialize_errc expected_err;

    error_sample(
        std::string&& name,
        deserialization_buffer&& from,
        metadata meta,
        deserialize_errc expected_err = deserialize_errc::protocol_value_error
    )
        : name(std::move(name)), from(std::move(from)), meta(std::move(meta)), expected_err(expected_err)
    {
    }
};

void add_int_samples(column_type type, std::size_t num_bytes, std::vector<error_sample>& output)
{
    output.emplace_back(error_sample(
        "signed_not_enough_space",
        deserialization_buffer(num_bytes, 0x0a),
        create_meta(type),
        deserialize_errc::incomplete_message
    ));
    output.emplace_back(error_sample(
        "unsigned_not_enough_space",
        deserialization_buffer(num_bytes, 0x0a),
        meta_builder().type(type).unsigned_flag(true).build(),
        deserialize_errc::incomplete_message
    ));
}

void add_bit_samples(std::vector<error_sample>& output)
{
    auto meta = meta_builder().type(column_type::bit).unsigned_flag(true).build();

    output.emplace_back(error_sample(
        "bit_error_deserializing_string_view",
        {0x01},
        meta,
        deserialize_errc::incomplete_message
    ));
    output.emplace_back(error_sample("bit_string_view_too_short", {0x00}, meta));
    output.emplace_back(error_sample(
        "bit_string_view_too_long",
        {0x09, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09},
        meta
    ));
}

void add_float_samples(std::vector<error_sample>& output)
{
    auto meta = create_meta(column_type::float_);

    output.push_back(
        error_sample("not_enough_space", {0x01, 0x02, 0x03}, meta, deserialize_errc::incomplete_message)
    );
    output.push_back(error_sample("inf", {0x00, 0x00, 0x80, 0x7f}, meta));
    output.push_back(error_sample("minus_inf", {0x00, 0x00, 0x80, 0xff}, meta));
    output.push_back(error_sample("nan", {0xff, 0xff, 0xff, 0x7f}, meta));
    output.push_back(error_sample("minus_nan", {0xff, 0xff, 0xff, 0xff}, meta));
}

void add_double_samples(std::vector<error_sample>& output)
{
    auto meta = create_meta(column_type::double_);

    output.push_back(error_sample(
        "not_enough_space",
        {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07},
        meta,
        deserialize_errc::incomplete_message
    ));
    output.push_back(error_sample("inf", {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x7f}, meta));
    output.push_back(error_sample("minus_inf", {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xff}, meta));
    output.push_back(error_sample("nan", {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f}, meta));
    output.push_back(error_sample("minus_nan", {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff}, meta));
}

// Based on correct, regular date {0x04, 0xda, 0x07, 0x03, 0x1c}
void add_date_samples(std::vector<error_sample>& output)
{
    auto meta = create_meta(column_type::date);

    output.push_back(error_sample("empty", {}, meta, deserialize_errc::incomplete_message));
    output.push_back(error_sample("incomplete_year", {0x04, 0xda}, meta, deserialize_errc::incomplete_message)
    );
    output.push_back(
        error_sample("no_month_day", {0x04, 0xda, 0x07}, meta, deserialize_errc::incomplete_message)
    );
    output.push_back(
        error_sample("no_day", {0x04, 0xda, 0x07, 0x03}, meta, deserialize_errc::incomplete_message)
    );
    output.push_back(error_sample(
        "invalid_year",
        {0x04, 0x10, 0x27, 0x03, 0x1c},  // year 10000
        meta
    ));
    output.push_back(error_sample("invalid_year_max", {0x04, 0xff, 0xff, 0x03, 0x1c}, meta));
    output.push_back(error_sample("invalid_month", {0x04, 0xda, 0x07, 13, 0x1c}, meta));
    output.push_back(error_sample("invalid_month_max", {0x04, 0xda, 0x07, 0xff, 0x1c}, meta));
    output.push_back(error_sample("invalid_day", {0x04, 0xda, 0x07, 0x03, 32}, meta));
    output.push_back(error_sample("invalid_day_max", {0x04, 0xda, 0x07, 0x03, 0xff}, meta));
    output.push_back(error_sample("protocol_max", {0xff, 0xff, 0xff, 0xff, 0xff}, meta));
}

// Based on correct datetime {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e,
// 0x00}
void add_datetime_samples(column_type type, std::vector<error_sample>& output)
{
    auto meta = create_meta(type);

    output.push_back(error_sample("empty", {}, meta, deserialize_errc::incomplete_message));
    output.push_back(
        error_sample("incomplete_date", {0x04, 0xda, 0x07, 0x01}, meta, deserialize_errc::incomplete_message)
    );
    output.push_back(error_sample(
        "no_hours_mins_secs",
        {0x07, 0xda, 0x07, 0x01, 0x01},
        meta,
        deserialize_errc::incomplete_message
    ));
    output.push_back(error_sample(
        "no_mins_secs",
        {0x07, 0xda, 0x07, 0x01, 0x01, 0x17},
        meta,
        deserialize_errc::incomplete_message
    ));
    output.push_back(error_sample(
        "no_secs",
        {0x07, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01},
        meta,
        deserialize_errc::incomplete_message
    ));
    output.push_back(error_sample(
        "incomplete_micros",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e},
        meta,
        deserialize_errc::incomplete_message
    ));
    output.push_back(error_sample("invalid_year_d", {0x04, 0x10, 0x27, 0x01, 0x01}, meta));  // year 10000
    output.push_back(error_sample("invalid_year_hms", {0x07, 0x10, 0x27, 0x01, 0x01, 0x17, 0x01, 0x3b}, meta)
    );
    output.push_back(error_sample(
        "invalid_year_hmsu",
        {0x0b, 0x10, 0x27, 0x01, 0x01, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        meta
    ));
    output.push_back(error_sample(
        "invalid_year_max_hmsu",
        {0x0b, 0xff, 0xff, 0x01, 0x01, 0x17, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        meta
    ));
    output.push_back(error_sample("invalid_hour_hms", {0x07, 0xda, 0x07, 0x01, 0x01, 24, 0x01, 0x3b}, meta));
    output.push_back(error_sample(
        "invalid_hour_hmsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 24, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        meta
    ));
    output.push_back(error_sample(
        "invalid_hour_max_hmsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0xff, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        meta
    ));
    output.push_back(error_sample("invalid_min_hms", {0x07, 0xda, 0x07, 0x01, 0x01, 0x17, 60, 0x3b}, meta));
    output.push_back(error_sample(
        "invalid_min_hmsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 60, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        meta
    ));
    output.push_back(error_sample(
        "invalid_min_max_hmsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0xff, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        meta
    ));
    output.push_back(error_sample("invalid_sec_hms", {0x07, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 60}, meta));
    output.push_back(error_sample(
        "invalid_sec_hmsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 60, 0x56, 0xc3, 0x0e, 0x00},
        meta
    ));
    output.push_back(error_sample(
        "invalid_sec_max_hmsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 0xff, 0x56, 0xc3, 0x0e, 0x00},
        meta
    ));
    output.push_back(error_sample(
        "invalid_micro_hmsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 0x3b, 0x40, 0x42, 0xf4, 0x00},
        meta
    ));  // 1M
    output.push_back(error_sample(
        "invalid_micro_max_hmsu",
        {0x0b, 0xda, 0x07, 0x01, 0x01, 0x17, 0x01, 0x3b, 0xff, 0xff, 0xff, 0xff},
        meta
    ));
    output.push_back(error_sample(
        "invalid_hour_invalid_date",
        {0x0b, 0x00, 0x00, 0x00, 0x00, 0xff, 0x01, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        meta
    ));
    output.push_back(error_sample(
        "invalid_min_invalid_date",
        {0x0b, 0x00, 0x00, 0x00, 0x00, 0x17, 0xff, 0x3b, 0x56, 0xc3, 0x0e, 0x00},
        meta
    ));
    output.push_back(error_sample(
        "invalid_sec_invalid_date",
        {0x0b, 0x00, 0x00, 0x00, 0x00, 0x17, 0x01, 0xff, 0x56, 0xc3, 0x0e, 0x00},
        meta
    ));
    output.push_back(error_sample(
        "invalid_micro_invalid_date",
        {0x0b, 0x00, 0x00, 0x00, 0x00, 0x17, 0x01, 0x3b, 0xff, 0xff, 0xff, 0xff},
        meta
    ));
    output.push_back(error_sample(
        "protocol_max",
        {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff},
        meta
    ));
}

void add_time_samples(std::vector<error_sample>& output)
{
    auto meta = create_meta(column_type::time);

    output.push_back(error_sample("empty", {}, meta, deserialize_errc::incomplete_message));
    output.push_back(
        error_sample("no_sign_days_hours_mins_secs", {0x08}, meta, deserialize_errc::incomplete_message)
    );
    output.push_back(
        error_sample("no_days_hours_mins_secs", {0x08, 0x01}, meta, deserialize_errc::incomplete_message)
    );
    output.push_back(error_sample(
        "no_hours_mins_secs",
        {0x08, 0x01, 0x22, 0x00, 0x00, 0x00},
        meta,
        deserialize_errc::incomplete_message
    ));
    output.push_back(error_sample(
        "no_mins_secs",
        {0x08, 0x01, 0x22, 0x00, 0x00, 0x00, 0x16},
        meta,
        deserialize_errc::incomplete_message
    ));
    output.push_back(error_sample(
        "no_secs",
        {0x08, 0x01, 0x22, 0x00, 0x00, 0x00, 0x16, 0x3b},
        meta,
        deserialize_errc::incomplete_message
    ));
    output.push_back(error_sample(
        "no_micros",
        {0x0c, 0x01, 0x22, 0x00, 0x00, 0x00, 0x16, 0x3b, 0x3a},
        meta,
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
        {"invalid_micros_max", {0x0c, 0x01, 0x22, 0x00, 0x00, 0x00, 0x16, 0x3b, 0x3a, 0xff, 0xff, 0xff, 0xff}
        },
    };

    for (auto& c : out_of_range_cases)
    {
        // Positive
        c.second[1] = 0x00;
        output.emplace_back(c.first + std::string("_positive"), deserialization_buffer(c.second), meta);

        // Negative
        c.second[1] = 0x01;
        output.emplace_back(c.first + std::string("_negative"), deserialization_buffer(c.second), meta);
    }
}

std::vector<error_sample> make_all_samples()
{
    std::vector<error_sample> res;
    add_int_samples(column_type::tinyint, 0, res);
    add_int_samples(column_type::smallint, 1, res);
    add_int_samples(column_type::mediumint, 3, res);
    add_int_samples(column_type::int_, 3, res);
    add_int_samples(column_type::bigint, 7, res);
    add_int_samples(column_type::year, 1, res);
    add_bit_samples(res);
    add_float_samples(res);
    add_double_samples(res);
    add_date_samples(res);
    add_datetime_samples(column_type::datetime, res);
    add_datetime_samples(column_type::timestamp, res);
    add_time_samples(res);
    return res;
}

BOOST_AUTO_TEST_CASE(error)
{
    for (const auto& tc : make_all_samples())
    {
        BOOST_TEST_CONTEXT("type=" << tc.meta.type() << ", name=" << tc.name)
        {
            detail::deserialization_context ctx(tc.from);

            field_view actual_value;
            auto err = deserialize_binary_field(ctx, tc.meta, actual_value);

            BOOST_TEST(err == tc.expected_err);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

}  // namespace
