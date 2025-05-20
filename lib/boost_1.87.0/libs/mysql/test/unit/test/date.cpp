//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/date.hpp>
#include <boost/mysql/days.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/detail/access.hpp>

#include <boost/core/span.hpp>
#include <boost/test/unit_test.hpp>

#include <chrono>
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

BOOST_AUTO_TEST_CASE(ctor_from_time_point_valid)
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

#ifdef BOOST_MYSQL_HAS_LOCAL_TIME
            date d2(std::chrono::local_days(std::chrono::days(tc.days_since_epoch)));
            BOOST_TEST(d2.valid());
            BOOST_TEST(d2.year() == tc.year);
            BOOST_TEST(d2.month() == tc.month);
            BOOST_TEST(d2.day() == tc.day);
#endif
        }
    }
}

BOOST_AUTO_TEST_CASE(ctor_from_time_point_invalid)
{
    BOOST_CHECK_THROW(date(date::time_point(days(2932897))), std::out_of_range);
    BOOST_CHECK_THROW(date(date::time_point(days(-719529))), std::out_of_range);
    BOOST_CHECK_THROW(date(date::time_point(days(INT_MAX))), std::out_of_range);
    BOOST_CHECK_THROW(date(date::time_point(days(INT_MIN))), std::out_of_range);
}

#ifdef BOOST_MYSQL_HAS_LOCAL_TIME
BOOST_AUTO_TEST_CASE(ctor_from_local_days_invalid)
{
    BOOST_CHECK_THROW(date(std::chrono::local_days(std::chrono::days(2932897))), std::out_of_range);
    BOOST_CHECK_THROW(date(std::chrono::local_days(std::chrono::days(-719529))), std::out_of_range);
    BOOST_CHECK_THROW(date(std::chrono::local_days(std::chrono::days(INT_MAX))), std::out_of_range);
    BOOST_CHECK_THROW(date(std::chrono::local_days(std::chrono::days(INT_MIN))), std::out_of_range);
}
#endif

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
    BOOST_TEST(date(2024, 5, 17).get_time_point().time_since_epoch().count() == 19860);
    BOOST_TEST(date(0, 1, 1).get_time_point().time_since_epoch().count() == -719528);
}

BOOST_AUTO_TEST_CASE(as_time_point)
{
    BOOST_TEST(date(9999, 12, 31).as_time_point().time_since_epoch().count() == 2932896);
    BOOST_TEST(date(0, 1, 1).as_time_point().time_since_epoch().count() == -719528);
    BOOST_TEST(date(2024, 5, 17).as_time_point().time_since_epoch().count() == 19860);
    BOOST_CHECK_THROW(date().as_time_point(), std::invalid_argument);
    BOOST_CHECK_THROW(date(2019, 2, 29).as_time_point(), std::invalid_argument);
}

#ifdef BOOST_MYSQL_HAS_LOCAL_TIME
BOOST_AUTO_TEST_CASE(get_local_time_point)
{
    BOOST_TEST(date(9999, 12, 31).get_local_time_point().time_since_epoch().count() == 2932896);
    BOOST_TEST(date(2024, 5, 17).as_local_time_point().time_since_epoch().count() == 19860);
    BOOST_TEST(date(0, 1, 1).get_local_time_point().time_since_epoch().count() == -719528);
}

BOOST_AUTO_TEST_CASE(as_local_time_point)
{
    BOOST_TEST(date(9999, 12, 31).as_local_time_point().time_since_epoch().count() == 2932896);
    BOOST_TEST(date(0, 1, 1).as_local_time_point().time_since_epoch().count() == -719528);
    BOOST_TEST(date(2024, 5, 17).as_local_time_point().time_since_epoch().count() == 19860);
    BOOST_CHECK_THROW(date().as_local_time_point(), std::invalid_argument);
    BOOST_CHECK_THROW(date(2019, 2, 29).as_local_time_point(), std::invalid_argument);
}
#endif

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

// operator stream is implemented in terms of to_string (see dt_to_string.hpp)
BOOST_AUTO_TEST_CASE(operator_stream)
{
    BOOST_TEST(stringize(date(2022, 1, 3)) == "2022-01-03");
    BOOST_TEST(stringize(date(2023, 12, 31)) == "2023-12-31");
    BOOST_TEST(stringize(date(0, 0, 0)) == "0000-00-00");
    BOOST_TEST(stringize(date(0xffff, 0xff, 0xff)) == "65535-255-255");
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
