//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/blob.hpp>
#include <boost/mysql/blob_view.hpp>
#include <boost/mysql/date.hpp>
#include <boost/mysql/datetime.hpp>
#include <boost/mysql/field.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/test/unit_test.hpp>

#include <cstdint>
#include <sstream>
#include <vector>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/create_basic.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

// Don't attempt to print std::chrono values
BOOST_TEST_DONT_PRINT_LOG_VALUE(boost::mysql::time)

BOOST_AUTO_TEST_SUITE(test_field_view)

BOOST_AUTO_TEST_SUITE(constructors)
BOOST_AUTO_TEST_CASE(default_constructor)
{
    field_view v;
    BOOST_TEST(v.is_null());
}

BOOST_AUTO_TEST_CASE(copy)
{
    field_view v(32);
    field_view v2(v);
    BOOST_TEST(v2.as_int64() == 32);
}

BOOST_AUTO_TEST_CASE(move)
{
    field_view v(field_view(32));
    BOOST_TEST(v.as_int64() == 32);
}

BOOST_AUTO_TEST_CASE(from_nullptr)
{
    field_view v(nullptr);
    BOOST_TEST(v.is_null());
}

BOOST_AUTO_TEST_CASE(from_u8)
{
    field_view v(std::uint8_t(0xfe));
    BOOST_TEST(v.as_uint64() == 0xfeu);
}

BOOST_AUTO_TEST_CASE(from_u16)
{
    field_view v(std::uint16_t(0xfefe));
    BOOST_TEST(v.as_uint64() == 0xfefeu);
}

BOOST_AUTO_TEST_CASE(from_u32)
{
    field_view v(std::uint32_t(0xfefefefe));
    BOOST_TEST(v.as_uint64() == 0xfefefefeu);
}

BOOST_AUTO_TEST_CASE(from_u64)
{
    field_view v(std::uint64_t(0xfefefefefefefefe));
    BOOST_TEST(v.as_uint64() == 0xfefefefefefefefeu);
}

BOOST_AUTO_TEST_CASE(from_s8)
{
    field_view v(std::int8_t(-1));
    BOOST_TEST(v.as_int64() == -1);
}

BOOST_AUTO_TEST_CASE(from_s16)
{
    field_view v(std::int16_t(-1));
    BOOST_TEST(v.as_int64() == -1);
}

BOOST_AUTO_TEST_CASE(from_s32)
{
    field_view v(std::int32_t(-1));
    BOOST_TEST(v.as_int64() == -1);
}

BOOST_AUTO_TEST_CASE(from_s64)
{
    field_view v(std::int64_t(-1));
    BOOST_TEST(v.as_int64() == -1);
}

BOOST_AUTO_TEST_CASE(from_char_array)
{
    field_view v("test");
    BOOST_TEST(v.as_string() == "test");
}

BOOST_AUTO_TEST_CASE(from_c_str)
{
    const char* str = "test";
    field_view v(str);
    BOOST_TEST(v.as_string() == "test");
}

BOOST_AUTO_TEST_CASE(from_string_view)
{
    string_view sv("test123", 4);
    field_view v(sv);
    BOOST_TEST(v.as_string() == "test");
}

BOOST_AUTO_TEST_CASE(from_blob_view)
{
    std::uint8_t buff[] = {0x00, 0x01, 0x02};
    field_view v{blob_view(buff)};
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(v.as_blob(), buff);
}

BOOST_AUTO_TEST_CASE(from_float)
{
    field_view v(4.2f);
    BOOST_TEST(v.as_float() == 4.2f);
}

BOOST_AUTO_TEST_CASE(from_double)
{
    field_view v(4.2);
    BOOST_TEST(v.as_double() == 4.2);
}

BOOST_AUTO_TEST_CASE(from_date)
{
    date d(2022, 4, 1);
    field_view v(d);
    BOOST_TEST(v.as_date() == d);
}

BOOST_AUTO_TEST_CASE(from_datetime)
{
    datetime d(2022u, 4u, 1u, 21u);
    field_view v(d);
    BOOST_TEST(v.as_datetime() == d);
}

BOOST_AUTO_TEST_CASE(from_time)
{
    auto t = maket(20, 10, 1);
    field_view v(t);
    BOOST_TEST(v.as_time() == t);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(accesors)

// Owning fields, to create references to them in the table below
field f_null;
field f_int64(50);
field f_uint64(50u);
field f_string("long_test_string");
field f_blob(blob{0x00, 0x01, 0x02, 0x03});
field f_float(4.2f);
field f_double(5.0);
field f_date(date(2020u, 1u, 1u));
field f_datetime(datetime(2019u, 1u, 1u));
field f_time(maket(9, 1, 0));

// clang-format off
struct
{
    const char* name;
    field_view field;
    field_kind expected_kind;
    bool is_null, is_int64, is_uint64, is_string, is_blob, is_float, is_double, is_date, is_datetime, is_time;
} test_cases [] = {
    // name           field                                kind                  null,  i64    u64    str    blob   float  double date   dt     time 
    { "null",         field_view(),                        field_kind::null,     true,  false, false, false, false, false, false, false, false, false },
    { "int64",        field_view(42),                      field_kind::int64,    false, true,  false, false, false, false, false, false, false, false },
    { "uint64",       field_view(42u),                     field_kind::uint64,   false, false, true,  false, false, false, false, false, false, false },
    { "string",       field_view("test"),                  field_kind::string,   false, false, false, true,  false, false, false, false, false, false },
    { "blob",         field_view(makebv("\0\1\xff")),      field_kind::blob,     false, false, false, false, true,  false, false, false, false, false },
    { "float",        field_view(4.2f),                    field_kind::float_,   false, false, false, false, false, true,  false, false, false, false },
    { "double",       field_view(4.2),                     field_kind::double_,  false, false, false, false, false, false, true,  false, false, false },
    { "date",         field_view(date(2020u, 1u, 1u)),     field_kind::date,     false, false, false, false, false, false, false, true,  false, false },
    { "datetime",     field_view(datetime(2020u, 1u, 1u)), field_kind::datetime, false, false, false, false, false, false, false, false, true,  false },
    { "time",         field_view(maket(20u, 1u, 1u)),      field_kind::time,     false, false, false, false, false, false, false, false, false, true },
    { "ref_null",     field_view(f_null),                  field_kind::null,     true,  false, false, false, false, false, false, false, false, false },
    { "ref_int64",    field_view(f_int64),                 field_kind::int64,    false, true,  false, false, false, false, false, false, false, false },
    { "ref_uint64",   field_view(f_uint64),                field_kind::uint64,   false, false, true,  false, false, false, false, false, false, false },
    { "ref_string",   field_view(f_string),                field_kind::string,   false, false, false, true,  false, false, false, false, false, false },
    { "ref_blob",     field_view(f_blob),                  field_kind::blob,     false, false, false, false, true,  false, false, false, false, false },
    { "ref_float",    field_view(f_float),                 field_kind::float_,   false, false, false, false, false, true,  false, false, false, false },
    { "ref_double",   field_view(f_double),                field_kind::double_,  false, false, false, false, false, false, true,  false, false, false },
    { "ref_date",     field_view(f_date),                  field_kind::date,     false, false, false, false, false, false, false, true,  false, false },
    { "ref_datetime", field_view(f_datetime),              field_kind::datetime, false, false, false, false, false, false, false, false, true,  false },
    { "ref_time",     field_view(f_time),                  field_kind::time,     false, false, false, false, false, false, false, false, false, true },
};
// clang-format on

BOOST_AUTO_TEST_CASE(kind)
{
    for (const auto& tc : test_cases)
    {
        BOOST_TEST(tc.field.kind() == tc.expected_kind, tc.name);
    }
}

BOOST_AUTO_TEST_CASE(is)
{
    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            BOOST_TEST(tc.field.is_null() == tc.is_null);
            BOOST_TEST(tc.field.is_int64() == tc.is_int64);
            BOOST_TEST(tc.field.is_uint64() == tc.is_uint64);
            BOOST_TEST(tc.field.is_string() == tc.is_string);
            BOOST_TEST(tc.field.is_blob() == tc.is_blob);
            BOOST_TEST(tc.field.is_float() == tc.is_float);
            BOOST_TEST(tc.field.is_double() == tc.is_double);
            BOOST_TEST(tc.field.is_date() == tc.is_date);
            BOOST_TEST(tc.field.is_datetime() == tc.is_datetime);
            BOOST_TEST(tc.field.is_time() == tc.is_time);
        }
    }
}

BOOST_AUTO_TEST_CASE(as_exceptions)
{
    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            if (tc.is_int64)
                BOOST_CHECK_NO_THROW(tc.field.as_int64());
            else
                BOOST_CHECK_THROW(tc.field.as_int64(), boost::mysql::bad_field_access);

            if (tc.is_uint64)
                BOOST_CHECK_NO_THROW(tc.field.as_uint64());
            else
                BOOST_CHECK_THROW(tc.field.as_uint64(), boost::mysql::bad_field_access);

            if (tc.is_string)
                BOOST_CHECK_NO_THROW(tc.field.as_string());
            else
                BOOST_CHECK_THROW(tc.field.as_string(), boost::mysql::bad_field_access);

            if (tc.is_blob)
                BOOST_CHECK_NO_THROW(tc.field.as_blob());
            else
                BOOST_CHECK_THROW(tc.field.as_blob(), boost::mysql::bad_field_access);

            if (tc.is_float)
                BOOST_CHECK_NO_THROW(tc.field.as_float());
            else
                BOOST_CHECK_THROW(tc.field.as_float(), boost::mysql::bad_field_access);

            if (tc.is_double)
                BOOST_CHECK_NO_THROW(tc.field.as_double());
            else
                BOOST_CHECK_THROW(tc.field.as_double(), boost::mysql::bad_field_access);

            if (tc.is_date)
                BOOST_CHECK_NO_THROW(tc.field.as_date());
            else
                BOOST_CHECK_THROW(tc.field.as_date(), boost::mysql::bad_field_access);

            if (tc.is_datetime)
                BOOST_CHECK_NO_THROW(tc.field.as_datetime());
            else
                BOOST_CHECK_THROW(tc.field.as_datetime(), boost::mysql::bad_field_access);

            if (tc.is_time)
                BOOST_CHECK_NO_THROW(tc.field.as_time());
            else
                BOOST_CHECK_THROW(tc.field.as_time(), boost::mysql::bad_field_access);
        }
    }
}

// Success cases (the type matches the called function)
BOOST_AUTO_TEST_CASE(int64)
{
    field_view f(-1);
    BOOST_TEST(f.as_int64() == -1);
    BOOST_TEST(f.get_int64() == -1);
}

BOOST_AUTO_TEST_CASE(uint64)
{
    field_view f(42u);
    BOOST_TEST(f.as_uint64() == 42u);
    BOOST_TEST(f.get_uint64() == 42u);
}

BOOST_AUTO_TEST_CASE(string)
{
    field_view f("test");
    BOOST_TEST(f.as_string() == "test");
    BOOST_TEST(f.get_string() == "test");
}

BOOST_AUTO_TEST_CASE(blob_)
{
    std::uint8_t buff[] = {0x00, 0x0f, 0x01};
    field_view f{blob_view(buff)};
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(f.as_blob(), buff);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(f.get_blob(), buff);
}

BOOST_AUTO_TEST_CASE(float_)
{
    field_view f(4.2f);
    BOOST_TEST(f.as_float() == 4.2f);
    BOOST_TEST(f.get_float() == 4.2f);
}

BOOST_AUTO_TEST_CASE(double_)
{
    field_view f(4.2);
    BOOST_TEST(f.as_double() == 4.2);
    BOOST_TEST(f.get_double() == 4.2);
}

BOOST_AUTO_TEST_CASE(date_)
{
    date d(2020u, 1u, 2u);
    field_view f(d);
    BOOST_TEST(f.as_date() == d);
    BOOST_TEST(f.get_date() == d);
}

BOOST_AUTO_TEST_CASE(datetime_)
{
    datetime dt(2020u, 1u, 2u);
    field_view f(dt);
    BOOST_TEST(f.as_datetime() == dt);
    BOOST_TEST(f.get_datetime() == dt);
}

BOOST_AUTO_TEST_CASE(time)
{
    auto t = maket(2020, 1, 2);
    field_view f(t);
    BOOST_TEST(f.as_time() == t);
    BOOST_TEST(f.get_time() == t);
}

BOOST_AUTO_TEST_CASE(ref_int64)
{
    field f(-1);
    field_view fv(f);
    BOOST_TEST(fv.as_int64() == -1);
    BOOST_TEST(fv.get_int64() == -1);
}

BOOST_AUTO_TEST_CASE(ref_uint64)
{
    field f(42u);
    field_view fv(f);
    BOOST_TEST(fv.as_uint64() == 42u);
    BOOST_TEST(fv.get_uint64() == 42u);
}

BOOST_AUTO_TEST_CASE(ref_string)
{
    field f("test");
    field_view fv(f);
    BOOST_TEST(fv.as_string() == "test");
    BOOST_TEST(fv.get_string() == "test");
}

BOOST_AUTO_TEST_CASE(ref_blob)
{
    std::uint8_t buff[] = {0x00, 0x01, 0x02};
    field f{blob(std::begin(buff), std::end(buff))};
    field_view fv(f);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(fv.as_blob(), buff);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(fv.get_blob(), buff);
}

BOOST_AUTO_TEST_CASE(ref_float)
{
    field f(4.2f);
    field_view fv(f);
    BOOST_TEST(fv.as_float() == 4.2f);
    BOOST_TEST(fv.get_float() == 4.2f);
}

BOOST_AUTO_TEST_CASE(ref_double)
{
    field f(4.2);
    field_view fv(f);
    BOOST_TEST(fv.as_double() == 4.2);
    BOOST_TEST(fv.get_double() == 4.2);
}

BOOST_AUTO_TEST_CASE(ref_date)
{
    date d(2020u, 1u, 2u);
    field f(d);
    field_view fv(f);
    BOOST_TEST(fv.as_date() == d);
    BOOST_TEST(fv.get_date() == d);
}

BOOST_AUTO_TEST_CASE(ref_datetime)
{
    datetime dt(2020u, 1u, 2u);
    field f(dt);
    field_view fv(f);
    BOOST_TEST(fv.as_datetime() == dt);
    BOOST_TEST(fv.get_datetime() == dt);
}

BOOST_AUTO_TEST_CASE(ref_time)
{
    auto t = maket(2020, 1, 2);
    field f(t);
    field_view fv(t);
    BOOST_TEST(fv.as_time() == t);
    BOOST_TEST(fv.get_time() == t);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_CASE(operator_equals)
{
    // clang-format off
    struct
    {
        const char* name;
        field_view f1;
        field_view f2;
        bool is_equal;
    } test_cases [] = {
        { "null_null", field_view(), field_view(), true, },
        { "null_int64", field_view(), field_view(-1), false },
        { "null_uint64", field_view(), field_view(42), false },
        { "null_string", field_view(), field_view("<NULL>"), false },
        { "null_blob", field_view(), field_view(blob_view()), false },
        { "null_float", field_view(), field_view(4.2f), false },
        { "null_double", field_view(), field_view(4.3), false },
        { "null_date", field_view(), field_view(date(2020u, 1u, 2u)), false },
        { "null_datetime", field_view(), field_view(datetime(2020u, 1u, 1u)), false },
        { "null_time", field_view(), field_view(maket(23, 1, 1)), false },

        { "int64_int64_same", field_view(42), field_view(42), true },
        { "int64_int64_different", field_view(42), field_view(-1), false },
        { "int64_uint64_same", field_view(42), field_view(42u), true },
        { "int64_uint64_different", field_view(42), field_view(43u), false },
        { "int64_uint64_zero", field_view(0), field_view(0u), true },
        { "int64_uint64_lt0", field_view(-1), field_view(42u), false },
        { "int64_uint64_gtmax", field_view(42), field_view(0xffffffffffffffffu), false },
        { "int64_uint64_lt0gtmax", field_view(-1), field_view(0xffffffffffffffffu), false },
        { "int64_string", field_view(42), field_view("42"), false },
        { "int64_blob", field_view(42), field_view(makebv("42")), false },
        { "int64_float", field_view(42), field_view(42.0f), false },
        { "int64_double", field_view(42), field_view(42.0), false },
        { "int64_date", field_view(42), field_view(date(2020u, 1u, 1u)), false },
        { "int64_datetime", field_view(42), field_view(datetime(2020u, 1u, 1u)), false },
        { "int64_time", field_view(42), field_view(maket(20, 1, 1)), false },

        { "uint64_uint64_same", field_view(0xffffffffffffffffu), field_view(0xffffffffffffffffu), true },
        { "uint64_uint64_different", field_view(42u), field_view(31u), false },
        { "uint64_string", field_view(42u), field_view("42"), false },
        { "uint64_blob", field_view(42u), field_view(makebv("42")), false },
        { "uint64_float", field_view(42u), field_view(42.0f), false },
        { "uint64_double", field_view(42u), field_view(42.0), false },
        { "uint64_date", field_view(42u), field_view(date(2020u, 1u, 1u)), false },
        { "uint64_datetime", field_view(42u), field_view(datetime(2020u, 1u, 1u)), false },
        { "uint64_time", field_view(42u), field_view(maket(20, 1, 1)), false },

        { "string_string_same", field_view("test"), field_view("test"), true },
        { "string_string_different", field_view("test"), field_view("test2"), false },
        { "string_blob", field_view("test"), field_view(makebv("test")), false },
        { "string_float", field_view("4.2"), field_view(4.2f), false },
        { "string_double", field_view("4.2"), field_view(4.2), false },
        { "string_date", field_view("2020-01-01"), field_view(date(2020u, 1u, 1u)), false },
        { "string_datetime", field_view("test"), field_view(datetime(2020u, 1u, 1u)), false },
        { "string_time", field_view("test"), field_view(maket(8, 1, 1)), false },
        
        { "blob_blob_same", field_view(makebv("\0test")), field_view(makebv("\0test")), true },
        { "blob_blob_different", field_view(makebv("\0test")), field_view(makebv("\0test2")), false },
        { "blob_blob_different_same_size", field_view(makebv("\0test")), field_view(makebv("\1test")), false },
        { "blob_blob_same_empty", field_view(blob_view()), field_view(blob_view()), true },
        { "blob_blob_different_empty", field_view(blob_view()), field_view(makebv("\0test")), false },
        { "blob_float", field_view(makebv("\x40\x86\x66\x66")), field_view(4.2f), false },
        { "blob_double", field_view(makebv("4.2")), field_view(4.2), false },
        { "blob_date", field_view(makebv("2020-01-01")), field_view(date(2020u, 1u, 1u)), false },
        { "blob_datetime", field_view(makebv("test")), field_view(datetime(2020u, 1u, 1u)), false },
        { "blob_time", field_view(makebv("test")), field_view(maket(8, 1, 1)), false },

        { "float_float_same", field_view(4.2f), field_view(4.2f), true },
        { "float_float_different", field_view(4.2f), field_view(0.0f), false },
        { "float_double", field_view(4.2f), field_view(4.2), false },
        { "float_date", field_view(4.2f), field_view(date(2020u, 1u, 2u)), false },
        { "float_datetime", field_view(4.2f), field_view(datetime(2020u, 1u, 2u)), false },
        { "float_time", field_view(4.2f), field_view(maket(20, 1, 2)), false },

        { "double_double_same", field_view(4.2), field_view(4.2), true },
        { "double_double_different", field_view(4.2), field_view(-1.0), false },
        { "double_date", field_view(4.2), field_view(date(2020u, 1u, 1u)), false },
        { "double_datetime", field_view(4.2), field_view(datetime(2020u, 1u, 1u)), false },
        { "double_time", field_view(4.2), field_view(maket(9, 1, 1)), false },

        { "date_date_same", field_view(date(2020u, 1u, 1u)), field_view(date(2020u, 1u, 1u)), true },
        { "date_date_different", field_view(date(2020u, 1u, 1u)), field_view(date(2019u, 1u, 1u)), false },
        { "date_datetime", field_view(date(2020u, 1u, 1u)), field_view(datetime(2020u, 1u, 1u)), false },
        { "date_time", field_view(date(2020u, 1u, 1u)), field_view(maket(9, 1, 1)), false },

        { "datetime_datetime_same", field_view(datetime(2020u, 1u, 1u, 10u)), field_view(datetime(2020u, 1u, 1u, 10u)), true },
        { "datetime_datetime_different", field_view(datetime(2020u, 1u, 1u, 10u)), field_view(datetime(2020u, 1u, 1u, 9u)), false },
        { "datetime_time", field_view(datetime(2020u, 1u, 1u)), field_view(maket(20, 1, 1)), false },

        { "time_time_same", field_view(maket(20, 1, 1)), field_view(maket(20, 1, 1)), true },
        { "time_time_different", field_view(maket(20, 1, 1)), field_view(maket(20, 1, 1, 10)), false },
    };
    // clang-format on

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // We compare regular field_view's and field_view's holding
            // pointers to fields, using the same cases to reduce duplication
            field owning_1(tc.f1);
            field owning_2(tc.f2);

            field_view f1 = tc.f1;
            field_view f2 = tc.f2;

            field_view fref1(owning_1);
            field_view fref2(owning_2);

            if (tc.is_equal)
            {
                BOOST_TEST(f1 == f2);
                BOOST_TEST(f2 == f1);
                BOOST_TEST(fref1 == fref2);
                BOOST_TEST(fref2 == fref1);
                BOOST_TEST(f1 == fref2);
                BOOST_TEST(fref2 == f1);

                BOOST_TEST(!(f1 != f2));
            }
            else
            {
                BOOST_TEST(!(f1 == f2));
                BOOST_TEST(!(f2 == f1));
                BOOST_TEST(!(fref1 == fref2));
                BOOST_TEST(!(fref2 == fref1));
                BOOST_TEST(!(f1 == fref2));
                BOOST_TEST(!(fref2 == f1));

                BOOST_TEST(tc.f1 != tc.f2);
                BOOST_TEST(tc.f2 != tc.f1);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(operator_equals_self_compare)
{
    // clang-format off
    struct
    {
        const char* name;
        field_view f;
    } test_cases [] = {
        { "null", field_view() },
        { "int64", field_view(40), },
        { "uint64", field_view(42u) },
        { "string", field_view("test") },
        { "blob", field_view(makebv("\0\0ab\1")) },
        { "float", field_view(4.2f) },
        { "double", field_view(5.0) },
        { "date", field_view(date(2020u, 1u, 1u)) },
        { "datetime", field_view(datetime(2020u, 1u, 1u)) },
        { "time", field_view(maket(8, 1, 1)) }
    };
    // clang-format on

    for (const auto& tc : test_cases)
    {
        // Regular field_view
        BOOST_TEST(tc.f == tc.f);

        // Reference to an owning field
        field owning_field(tc.f);
        field_view fref(owning_field);
        BOOST_TEST(fref == fref);
    }
}

// operator<<
struct stream_sample
{
    std::string name;
    field_view input;
    std::string expected;

    template <class T>
    stream_sample(std::string&& name, const T& input, std::string&& expected)
        : name(std::move(name)), input(input), expected(std::move(expected))
    {
    }
};

BOOST_AUTO_TEST_CASE(operator_stream)
{
    struct owning_fields_t
    {
        field f_null{};
        field f_int64{-1};
        field f_uint64{50u};
        field f_string{"long_test_string"};
        field f_blob{
            blob{0x00, 0x32, 0x01}
        };
        field f_float{4.2f};
        field f_double{5.1};
        field f_date{date(2020u, 1u, 1u)};
        field f_datetime{datetime(2019u, 1u, 1u, 21u, 19u, 1u, 9u)};
        field f_time{maket(9, 1, 0, 210)};
    };
    static owning_fields_t owning_fields;

    struct test_case_t
    {
        string_view name;
        field_view input;
        string_view expected;
    } test_cases[] = {
        // Field views holding values
        {"null", field_view(), "<NULL>"},
        {"i64_positive", field_view(std::int64_t(42)), "42"},
        {"i64_negative", field_view(std::int64_t(-90)), "-90"},
        {"i64_zero", field_view(std::int64_t(0)), "0"},
        {"u64_positive", field_view(std::uint64_t(42)), "42"},
        {"u64_zero", field_view(std::uint64_t(0)), "0"},
        {"string_view", field_view("a_string"), "a_string"},
        {"blob_empty", field_view(blob_view()), "{}"},
        {"blob_one_elm", field_view(makebv("\4")), "{ 0x04 }"},
        {"blob_two_elms", field_view(makebv("\4\5")), "{ 0x04, 0x05 }"},
        {"blob_several_elms", field_view(makebv("\0\x0a\x2f\xff")), "{ 0x00, 0x0a, 0x2f, 0xff }"},
        {"blob_padding",
         field_view(makebv("\x0e\x0f\x0a\x0b\x1f\x20\x7f\x80\xff")),
         "{ 0x0e, 0x0f, 0x0a, 0x0b, 0x1f, 0x20, 0x7f, 0x80, 0xff }"},
        {"float", field_view(2.43f), "2.43"},
        {"double", field_view(8.12), "8.12"},
        {"date", field_view(date(2020, 1, 19)), "2020-01-19"},
        {"datetime", field_view(datetime(2020, 1, 19, 11, 30, 21, 98765)), "2020-01-19 11:30:21.098765"},
        {"time", field_view(-maket(12, 1, 5, 345)), "-12:01:05.000345"},

        // Field views holding pointers to fields
        {"ref_null", field_view(owning_fields.f_null), "<NULL>"},
        {"ref_int64", field_view(owning_fields.f_int64), "-1"},
        {"ref_uint64", field_view(owning_fields.f_uint64), "50"},
        {"ref_string", field_view(owning_fields.f_string), "long_test_string"},
        {"ref_blob", field_view(owning_fields.f_blob), "{ 0x00, 0x32, 0x01 }"},
        {"ref_float", field_view(owning_fields.f_float), "4.2"},
        {"ref_double", field_view(owning_fields.f_double), "5.1"},
        {"ref_date", field_view(owning_fields.f_date), "2020-01-01"},
        {"ref_datetime", field_view(owning_fields.f_datetime), "2019-01-01 21:19:01.000009"},
        {"ref_time", field_view(owning_fields.f_time), "09:01:00.000210"},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            std::ostringstream ss;
            ss << tc.input;
            BOOST_TEST(ss.str() == tc.expected);
        }
    }
}

// Make sure constxpr can actually be used in a constexpr context
// C++14+ only
#ifndef BOOST_NO_CXX14_CONSTEXPR

BOOST_AUTO_TEST_SUITE(constexpr_fns)

BOOST_AUTO_TEST_CASE(null)
{
    constexpr field_view v{};
    constexpr field_view v2(nullptr);
    static_assert(v.is_null(), "");
    static_assert(v == v2, "");
    static_assert(!(v != v2), "");
}

BOOST_AUTO_TEST_CASE(int64)
{
    constexpr field_view v(60);
    static_assert(v.is_int64(), "");
    static_assert(v.as_int64() == 60, "");
    static_assert(v.get_int64() == 60, "");
    static_assert(v == field_view(60), "");
}

BOOST_AUTO_TEST_CASE(uint64)
{
    constexpr field_view v(60u);
    static_assert(v.is_uint64(), "");
    static_assert(v.as_uint64() == 60u, "");
    static_assert(v.get_uint64() == 60u, "");
    static_assert(v == field_view(60u), "");
}

BOOST_AUTO_TEST_CASE(string)
{
    constexpr field_view v(makesv("test"));
    static_assert(v.is_string(), "");
    static_assert(v.as_string()[0] == 't', "");
    static_assert(v.get_string()[0] == 't', "");
    // string_view comparison uses char_traits::compare, which is not constexpr in C++14
}

BOOST_AUTO_TEST_CASE(blob)
{
    constexpr field_view v{blob_view()};
    static_assert(v.is_blob(), "");
    static_assert(v.as_blob().empty(), "");
    static_assert(v.get_blob().empty(), "");
    // blob comparison is not constexpr
}

BOOST_AUTO_TEST_CASE(float_)
{
    constexpr field_view v(4.2f);
    static_assert(v.is_float(), "");
    static_assert(v.as_float() == 4.2f, "");
    static_assert(v.get_float() == 4.2f, "");
    static_assert(v == field_view(4.2f), "");
}

BOOST_AUTO_TEST_CASE(double_)
{
    constexpr field_view v(4.2);
    static_assert(v.is_double(), "");
    static_assert(v.as_double() == 4.2, "");
    static_assert(v.get_double() == 4.2, "");
    static_assert(v == field_view(4.2), "");
}

BOOST_AUTO_TEST_CASE(date_)
{
    constexpr date d(2020, 1, 1);
    constexpr field_view v(d);
    static_assert(v.is_date(), "");
    static_assert(v.as_date() == d, "");
    static_assert(v.get_date() == d, "");
    static_assert(v == field_view(d), "");
}

BOOST_AUTO_TEST_CASE(datetime_)
{
    constexpr datetime d(2020u, 1u, 1u);
    constexpr field_view v(d);
    static_assert(v.is_datetime(), "");
    static_assert(v.as_datetime() == d, "");
    static_assert(v.get_datetime() == d, "");
    static_assert(v == field_view(d), "");
}

BOOST_AUTO_TEST_CASE(time)
{
    constexpr auto t = maket(8, 1, 1);
    constexpr field_view v(t);
    static_assert(v.is_time(), "");
    static_assert(v.as_time() == t, "");
    static_assert(v.get_time() == t, "");
    static_assert(v == field_view(t), "");
}

BOOST_AUTO_TEST_SUITE_END()

#endif  // BOOST_NO_CXX14_CONSTEXPR

BOOST_AUTO_TEST_SUITE_END()
