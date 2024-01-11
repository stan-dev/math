//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/detail/datetime.hpp>

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <string>

#include "test_common/stringize.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql::detail;

// These tests are very extensive in range. Making them parameterized
// proves very runtime expensive. We rather use plain loops. Using plain
// BOOST_TEST() also increases runtime way too much

BOOST_AUTO_TEST_SUITE(test_datetime_detail)

using test_name_t = std::array<char, 32>;

class test_state
{
    test_name_t test_name_{};
    std::vector<std::string> failed_assertions_;

public:
    const std::vector<std::string>& failures() const noexcept { return failed_assertions_; }
    void set_test_name(test_name_t v) noexcept { test_name_ = v; }

    template <class T1, class T2>
    void assert_equals(const T1& v1, const T2& v2, int line)
    {
        if (v1 != v2)
        {
            failed_assertions_.push_back(
                stringize(__FILE__, ":", line, " (context=", test_name_.data(), "): ", v1, " != ", v2)
            );
        }
    }

    void check()
    {
        for (const auto& fail : failed_assertions_)
        {
            BOOST_TEST(false, fail);
        }
    }
};

#define BOOST_MYSQL_ASSERT_EQ(st, v1, v2) st.assert_equals(v1, v2, __LINE__)

// Helpers
constexpr std::uint16_t leap_years[] = {
    1804, 1808, 1812, 1816, 1820, 1824, 1828, 1832, 1836, 1840, 1844, 1848, 1852, 1856, 1860, 1864, 1868,
    1872, 1876, 1880, 1884, 1888, 1892, 1896, 1904, 1908, 1912, 1916, 1920, 1924, 1928, 1932, 1936, 1940,
    1944, 1948, 1952, 1956, 1960, 1964, 1968, 1972, 1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008,
    2012, 2016, 2020, 2024, 2028, 2032, 2036, 2040, 2044, 2048, 2052, 2056, 2060, 2064, 2068, 2072, 2076,
    2080, 2084, 2088, 2092, 2096, 2104, 2108, 2112, 2116, 2120, 2124, 2128, 2132, 2136, 2140, 2144, 2148,
    2152, 2156, 2160, 2164, 2168, 2172, 2176, 2180, 2184, 2188, 2192, 2196, 2204,
};

bool is_leap_year(std::uint16_t y)
{
    return std::binary_search(std::begin(leap_years), std::end(leap_years), y);
}

std::uint8_t last_day_of_month(std::uint8_t month)  // doesn't take leap years into account
{
    constexpr std::uint8_t last_month_days[] = {31u, 28u, 31u, 30u, 31u, 30u, 31u, 31u, 30u, 31u, 30u, 31u};
    assert(month >= 1 && month <= 12);
    return last_month_days[month - 1];
}

std::array<char, 32> date_to_string(std::uint16_t year, std::uint8_t month, std::uint8_t day)
{
    std::array<char, 32> res{};
    snprintf(
        res.data(),
        res.size(),
        "%4u-%02u-%2u",
        static_cast<unsigned>(year),
        static_cast<unsigned>(month),
        static_cast<unsigned>(day)
    );
    return res;
}

BOOST_AUTO_TEST_SUITE(is_valid_)

// thorough coverage for 400 years
BOOST_AUTO_TEST_CASE(coverage)
{
    test_state st;

    for (std::uint16_t year = 1804; year <= 2204; ++year)
    {
        bool leap = is_leap(year);
        for (std::uint8_t month = 1; month <= 12; ++month)
        {
            std::uint8_t last_month_day = month == 2 && leap ? 29u : last_day_of_month(month);
            for (std::uint8_t day = 1; day <= 32; ++day)
            {
                st.set_test_name(date_to_string(year, month, day));

                BOOST_MYSQL_ASSERT_EQ(st, (is_valid(year, month, day)), (day <= last_month_day));
            }
        }
    }

    st.check();
}

// spotchecks for certain invalid dates
BOOST_AUTO_TEST_CASE(invalid_spotchecks)
{
    // year out of range of MySQL validity
    BOOST_TEST(!is_valid(10000u, 1u, 1u));
    BOOST_TEST(!is_valid(0xffffu, 1u, 1u));

    // month out of range
    BOOST_TEST(!is_valid(2010u, 13u, 1u));
    BOOST_TEST(!is_valid(2010u, 0u, 1u));
    BOOST_TEST(!is_valid(2010u, 0xffu, 1u));

    // day out of range
    BOOST_TEST(!is_valid(2019u, 2u, 29u));
    BOOST_TEST(!is_valid(2010u, 2u, 32u));
    BOOST_TEST(!is_valid(2010u, 2u, 0u));
    BOOST_TEST(!is_valid(2010u, 2u, 0xffu));

    // combinations
    BOOST_TEST(!is_valid(0u, 0u, 0u));
    BOOST_TEST(!is_valid(0xffffu, 0xffu, 0xffu));
    BOOST_TEST(!is_valid(2010u, 0u, 0u));
    BOOST_TEST(!is_valid(0xffffu, 42u, 0xffu));
}

// spotchecks for certain valid dates
BOOST_AUTO_TEST_CASE(valid_spotchecks)
{
    BOOST_TEST(is_valid(0u, 1u, 1u));
    BOOST_TEST(is_valid(2020u, 2u, 29u));
    BOOST_TEST(is_valid(9999u, 1u, 1u));
}

BOOST_AUTO_TEST_SUITE_END()

// ymd_to_days, days_to_ymd
// Helper function that actually performs the assertions for us
void ymd_years_test(test_state& st, std::uint16_t year, std::uint8_t month, std::uint8_t day, int num_days)
{
    st.set_test_name(date_to_string(year, month, day));

    BOOST_MYSQL_ASSERT_EQ(st, is_valid(year, month, day), true);
    BOOST_MYSQL_ASSERT_EQ(st, ymd_to_days(year, month, day), num_days);
    std::uint16_t output_year{};
    std::uint8_t output_month{}, output_day{};
    bool ok = days_to_ymd(num_days, output_year, output_month, output_day);

    BOOST_MYSQL_ASSERT_EQ(st, ok, true);
    BOOST_MYSQL_ASSERT_EQ(st, output_day, day);
    BOOST_MYSQL_ASSERT_EQ(st, output_month, month);
    BOOST_MYSQL_ASSERT_EQ(st, output_year, year);
}

BOOST_AUTO_TEST_CASE(ymd_to_days_days_to_ymd)
{
    test_state st;

    // Starting from 1970, going up
    int num_days = 0;

    for (int year = 1970; year <= 2204; ++year)
    {
        for (unsigned month = 1; month <= 12; ++month)
        {
            unsigned last_month_day = month == 2 && is_leap_year(year) ? 29u : last_day_of_month(month);
            for (unsigned day = 1; day <= last_month_day; ++day)
            {
                ymd_years_test(st, year, month, day, num_days++);
            }
        }
    }

    // Starting from 1970, going down
    num_days = -1;

    for (int year = 1969; year >= 1804; --year)
    {
        for (unsigned month = 12; month >= 1; --month)
        {
            unsigned last_month_day = month == 2 && is_leap_year(year) ? 29u : last_day_of_month(month);
            for (unsigned day = last_month_day; day >= 1; --day)
            {
                ymd_years_test(st, year, month, day, num_days--);
            }
        }
    }

    st.check();
}

BOOST_AUTO_TEST_CASE(ymd_to_days_spotcheck)
{
    BOOST_TEST(ymd_to_days(0, 1, 1) == -719528);
    BOOST_TEST(ymd_to_days(1970, 1, 1) == 0);
    BOOST_TEST(ymd_to_days(9999, 12, 31) == 2932896);
}

// Verify range checks work
BOOST_AUTO_TEST_CASE(days_to_ymd_limits)
{
    // Just in the lower limit
    std::uint16_t years;
    std::uint8_t month, day;
    bool ok = days_to_ymd(-719528, years, month, day);
    BOOST_TEST(ok);
    BOOST_TEST(years == 0);
    BOOST_TEST(month == 1);
    BOOST_TEST(day == 1);

    // Below lower limit
    BOOST_TEST(!days_to_ymd(-719529, years, month, day));
    BOOST_TEST(!days_to_ymd((std::numeric_limits<int>::min)(), years, month, day));

    // Just in the upper limit
    ok = days_to_ymd(2932896, years, month, day);
    BOOST_TEST(ok);
    BOOST_TEST(years == 9999);
    BOOST_TEST(month == 12);
    BOOST_TEST(day == 31);

    // Above the upper limit. 719468 is a magic number used within the algorithm that
    // was found to cause signed int overflow
    BOOST_TEST(!days_to_ymd(2932897, years, month, day));
    BOOST_TEST(!days_to_ymd((std::numeric_limits<int>::max)() - 719467, years, month, day));
    BOOST_TEST(!days_to_ymd((std::numeric_limits<int>::max)() - 719468, years, month, day));
    BOOST_TEST(!days_to_ymd((std::numeric_limits<int>::max)() - 719469, years, month, day));
    BOOST_TEST(!days_to_ymd((std::numeric_limits<int>::max)(), years, month, day));
}

BOOST_AUTO_TEST_SUITE_END()
