//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/date.hpp>
#include <boost/mysql/days.hpp>

#include <boost/test/unit_test.hpp>

#include <stdexcept>

#include "test_common/stringize.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

// Note: most of the algorithms have been thoroughly covered in
// detail/auxiliar/datetime.cpp, so just spotchecks here

BOOST_AUTO_TEST_SUITE(test_date)

BOOST_AUTO_TEST_CASE(default_constructor)
{
    date d;
    BOOST_TEST(d.year() == 0);
    BOOST_TEST(d.month() == 0);
    BOOST_TEST(d.day() == 0);
    BOOST_TEST(!d.valid());
}

BOOST_AUTO_TEST_SUITE(ctor_from_time_point)
BOOST_AUTO_TEST_CASE(valid)
{
    struct
    {
        int days_since_epoch;
        std::uint16_t year;
        std::uint8_t month;
        std::uint8_t day;
    } test_cases[] = {
        {2932896, 9999, 12, 31},
        {0,       1970, 1,  1 },
        {-719528, 0,    1,  1 },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.days_since_epoch)
        {
            date::time_point tp(days(tc.days_since_epoch));
            date d(tp);
            BOOST_TEST(d.valid());
            BOOST_TEST(d.year() == tc.year);
            BOOST_TEST(d.month() == tc.month);
            BOOST_TEST(d.day() == tc.day);
        }
    }
}

BOOST_AUTO_TEST_CASE(invalid)
{
    BOOST_CHECK_THROW(date(date::time_point(days(2932897))), std::out_of_range);
    BOOST_CHECK_THROW(date(date::time_point(days(-719529))), std::out_of_range);
    BOOST_CHECK_THROW(date(date::time_point(days(INT_MAX))), std::out_of_range);
    BOOST_CHECK_THROW(date(date::time_point(days(INT_MIN))), std::out_of_range);
}
BOOST_AUTO_TEST_SUITE_END()  // ctor_from_time_point

BOOST_AUTO_TEST_CASE(valid)
{
    BOOST_TEST(date(2020, 1, 1).valid());
    BOOST_TEST(date(2020, 2, 29).valid());
    BOOST_TEST(!date(2019, 2, 29).valid());
    BOOST_TEST(!date(0xffff, 0xff, 0xff).valid());
}

BOOST_AUTO_TEST_CASE(get_time_point)
{
    BOOST_TEST(date(9999, 12, 31).get_time_point().time_since_epoch().count() == 2932896);
    BOOST_TEST(date(0, 1, 1).get_time_point().time_since_epoch().count() == -719528);
}

BOOST_AUTO_TEST_CASE(as_time_point)
{
    BOOST_TEST(date(9999, 12, 31).as_time_point().time_since_epoch().count() == 2932896);
    BOOST_TEST(date(0, 1, 1).as_time_point().time_since_epoch().count() == -719528);
    BOOST_CHECK_THROW(date().as_time_point(), std::invalid_argument);
    BOOST_CHECK_THROW(date(2019, 2, 29).as_time_point(), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(operator_equals)
{
    struct
    {
        const char* name;
        date d1;
        date d2;
        bool equal;
    } test_cases[] = {
        {"equal",           date(2020,   2,    29),   date(2020,   2,    29),   true },
        {"equal_invalid",   date(0,      0,    0),    date(0,      0,    0),    true },
        {"equal_max",       date(0xffff, 0xff, 0xff), date(0xffff, 0xff, 0xff), true },
        {"different_year",  date(2020,   1,    31),   date(2019,   1,    31),   false},
        {"different_month", date(2020,   1,    10),   date(2020,   2,    10),   false},
        {"different_day",   date(2020,   1,    10),   date(2020,   1,    11),   false},
        {"all_different",   date(2020,   1,    1),    date(2021,   2,    2),    false},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            BOOST_TEST((tc.d1 == tc.d2) == tc.equal);
            BOOST_TEST((tc.d2 == tc.d1) == tc.equal);
            BOOST_TEST((tc.d1 != tc.d2) == !tc.equal);
            BOOST_TEST((tc.d2 != tc.d1) == !tc.equal);
        }
    }
}

BOOST_AUTO_TEST_CASE(operator_stream)
{
    struct
    {
        const char* name;
        std::uint16_t value;
        const char* repr;
    } year_values[] = {
        {"min",      0,      "0000" },
        {"onedig",   1,      "0001" },
        {"twodig",   98,     "0098" },
        {"threedig", 789,    "0789" },
        {"regular",  1999,   "1999" },
        {"fourdig",  9999,   "9999" },
        {"max",      0xffff, "65535"},
    };

    struct
    {
        const char* name;
        std::uint8_t value;
        const char* repr;
    } month_values[] = {
        {"zero", 0,    "00" },
        {"1dig", 2,    "02" },
        {"2dig", 12,   "12" },
        {"max",  0xff, "255"},
    };

    struct
    {
        const char* name;
        std::uint8_t value;
        const char* repr;
    } day_values[] = {
        {"zero", 0,    "00" },
        {"1dig", 1,    "01" },
        {"2dig", 31,   "31" },
        {"max",  0xff, "255"},
    };

    for (const auto& year : year_values)
    {
        for (const auto& month : month_values)
        {
            for (const auto& day : day_values)
            {
                BOOST_TEST_CONTEXT("year=" << year.name << ", month=" << month.name << ", day=" << day.name)
                {
                    auto expected = stringize(year.repr, '-', month.repr, '-', day.repr);
                    auto actual = stringize(date(year.value, month.value, day.value));
                    BOOST_TEST(actual == expected);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(now)
{
    auto d = date::now();
    BOOST_TEST_REQUIRE(d.valid());
    BOOST_TEST(d.year() > 2020u);
    BOOST_TEST(d.year() < 2100u);
}

// Make sure constxpr can actually be used in a constexpr context
BOOST_AUTO_TEST_CASE(constexpr_fns_cxx11)
{
    constexpr date d0{};
    static_assert(!d0.valid(), "");
    static_assert(d0.year() == 0, "");
    static_assert(d0.month() == 0, "");
    static_assert(d0.day() == 0, "");

    constexpr date d1(2020, 10, 1);
    static_assert(d1.valid(), "");
    static_assert(d1.year() == 2020, "");
    static_assert(d1.month() == 10, "");
    static_assert(d1.day() == 1, "");

    static_assert(d0 == date(), "");
    static_assert(d0 != d1, "");
}

#ifndef BOOST_NO_CXX14_CONSTEXPR
BOOST_AUTO_TEST_CASE(constexpr_fns_cxx14)
{
    constexpr date d0(date::time_point(days(2932896)));
    static_assert(d0 == date(9999, 12, 31), "");

    constexpr auto tp1 = d0.get_time_point();
    static_assert(tp1 == date::time_point(days(2932896)), "");

    constexpr auto tp2 = d0.as_time_point();
    static_assert(tp2 == tp1, "");
}
#endif

BOOST_AUTO_TEST_SUITE_END()
