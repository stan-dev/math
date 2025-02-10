//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/string_view.hpp>
#include <boost/mysql/time.hpp>

#include <boost/mysql/impl/internal/dt_to_string.hpp>

#include <boost/core/span.hpp>
#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <cstdint>
#include <random>

#include "test_common/create_basic.hpp"
#include "test_common/stringize.hpp"

// These tests use heap-allocated buffers for all cases because it helps asan detect overflows

using namespace boost::mysql;
using test::maket;

BOOST_AUTO_TEST_SUITE(test_dt_to_string)

// A wrapper around std::random to generate integers using a uniform distribution
template <class IntType>
class int_generator
{
    std::mt19937 mt_{std::random_device{}()};
    std::uniform_int_distribution<IntType> dist_;

public:
    int_generator(IntType a, IntType b) : dist_(a, b) {}
    IntType generate() { return dist_(mt_); }
};

// Represents a broken date/datetime/time component. Used in all coverage tests
template <class IntType>
struct component_value
{
    const char* name;
    IntType value;
    const char* repr;
};

// date
static std::string invoke_to_string(date d)
{
    std::string res(32, '\0');
    std::size_t sz = detail::date_to_string(d.year(), d.month(), d.day(), boost::span<char, 32>(&res[0], 32));
    res.resize(sz);
    return res;
}

BOOST_AUTO_TEST_CASE(date_coverage)
{
    // We use a dot product to cover common cases
    component_value<std::uint16_t> year_values[] = {
        {"min",      0,      "0000" },
        {"onedig",   1,      "0001" },
        {"twodig",   98,     "0098" },
        {"threedig", 789,    "0789" },
        {"regular",  1999,   "1999" },
        {"fourdig",  9999,   "9999" },
        {"max",      0xffff, "65535"},
    };

    component_value<std::uint8_t> month_values[] = {
        {"zero", 0,    "00" },
        {"1dig", 2,    "02" },
        {"2dig", 12,   "12" },
        {"max",  0xff, "255"},
    };

    component_value<std::uint8_t> day_values[] = {
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
                    // Expected value
                    auto expected = test::stringize(year.repr, '-', month.repr, '-', day.repr);

                    // Input value
                    date d(year.value, month.value, day.value);

                    // Call the function. Using a heap-allocated buffer helps detect overruns
                    auto actual = invoke_to_string(d);

                    // Check
                    BOOST_TEST(actual == expected);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(date_padding)
{
    // Double-check we correctly pad, regardless of the number
    // All dates below 9999-xx-xx should have 10 characters
    constexpr std::size_t expected_size = 10;

    // Day
    for (std::uint8_t day = 0u; day <= 31u; ++day)
    {
        BOOST_TEST_CONTEXT(day)
        {
            date d(2021, 1, day);
            BOOST_TEST(invoke_to_string(d).size() == expected_size);
        }
    }

    // Month
    for (std::uint8_t month = 0u; month <= 12u; ++month)
    {
        BOOST_TEST_CONTEXT(month)
        {
            date d(2021, month, 12);
            BOOST_TEST(invoke_to_string(d).size() == expected_size);
        }
    }

    // Year. Iterating all values is too costly, so we check some random ones
    int_generator<std::uint16_t> year_gen(std::uint16_t(0), std::uint16_t(9999));
    for (int i = 0; i < 30; ++i)
    {
        std::uint16_t year = year_gen.generate();
        BOOST_TEST_CONTEXT(year)
        {
            date d(year, 2, 12);
            BOOST_TEST(invoke_to_string(d).size() == expected_size);
        }
    }
}

// datetime
static std::string invoke_to_string(datetime d)
{
    std::string res(64, '\0');
    std::size_t sz = detail::datetime_to_string(
        d.year(),
        d.month(),
        d.day(),
        d.hour(),
        d.minute(),
        d.second(),
        d.microsecond(),
        boost::span<char, 64>(&res[0], 64)
    );
    res.resize(sz);
    return res;
}

BOOST_AUTO_TEST_CASE(datetime_coverage)
{
    // We use a dot product to cover common cases
    constexpr component_value<std::uint16_t> year_values[] = {
        {"min",       0,      "0000" },
        {"onedig",    1,      "0001" },
        {"twodig",    98,     "0098" },
        {"threedig",  789,    "0789" },
        {"regular",   1999,   "1999" },
        {"max_mysql", 9999,   "9999" },
        {"max",       0xffff, "65535"},
    };

    constexpr component_value<std::uint8_t> month_values[] = {
        {"zero",   0,    "00" },
        {"onedig", 2,    "02" },
        {"twodig", 12,   "12" },
        {"max",    0xff, "255"},
    };

    constexpr component_value<std::uint8_t> day_values[] = {
        {"zero",   0,    "00" },
        {"onedig", 1,    "01" },
        {"twodig", 31,   "31" },
        {"max",    0xff, "255"},
    };

    constexpr component_value<std::uint8_t> hours_values[] = {
        {"zero",   0,    "00" },
        {"onedig", 5,    "05" },
        {"twodig", 23,   "23" },
        {"max",    0xff, "255"},
    };

    constexpr component_value<std::uint8_t> mins_secs_values[] = {
        {"zero",   0,    "00" },
        {"onedig", 5,    "05" },
        {"twodig", 59,   "59" },
        {"max",    0xff, "255"},
    };

    constexpr component_value<std::uint32_t> micros_values[] = {
        {"zero",      0,          "000000"    },
        {"onedig",    5,          "000005"    },
        {"twodig",    50,         "000050"    },
        {"max_mysql", 999999,     "999999"    },
        {"max",       0xffffffff, "4294967295"},
    };

    // clang-format off
    for (const auto& year : year_values)
    {
    for (const auto& month : month_values)
    {
    for (const auto& day : day_values)
    {
    for (const auto& hours : hours_values)
    {
    for (const auto& mins : mins_secs_values)
    {
    for (const auto& secs : mins_secs_values)
    {
    for (const auto& micros : micros_values)
    {
        BOOST_TEST_CONTEXT(
            "year=" << year.name << ", month=" << month.name << "day=" << day.name <<
            "hour=" << hours.name << ", mins=" << mins.name << ", secs=" << secs.name <<
            "micros=" << micros.name
        )
        {
            // Expected value
            std::string expected = test::stringize(
                year.repr, '-', month.repr, '-', day.repr, ' ',
                hours.repr, ':', mins.repr, ':', secs.repr,
                '.', micros.repr
            );

            // Input value
            datetime dt(
                year.value,
                month.value,
                day.value,
                hours.value,
                mins.value,
                secs.value,
                micros.value
            );

            // Call the function
            auto actual = invoke_to_string(dt);

            // Check
            BOOST_TEST(actual == expected);
        }
    }
    }
    }
    }
    }
    }
    }
    // clang-format on
}

BOOST_AUTO_TEST_CASE(datetime_padding)
{
    // Double-check we correctly pad, regardless of the number
    // All datetimes below 9999-xx-xx xx:xx:xx.xxxxxx should have 26 characters
    constexpr std::size_t expected_size = 26;

    // Year. Iterating all values is too costly, so we check some random ones
    int_generator<std::uint16_t> year_gen(std::uint16_t(0), std::uint16_t(9999));
    for (int i = 0; i < 30; ++i)
    {
        std::uint16_t year = year_gen.generate();
        BOOST_TEST_CONTEXT(year)
        {
            datetime d(year, 2, 12);
            BOOST_TEST(invoke_to_string(d).size() == expected_size);
        }
    }

    // Month
    for (std::uint8_t month = 0u; month <= 31u; ++month)
    {
        BOOST_TEST_CONTEXT(month)
        {
            datetime d(2021, month, 12);
            BOOST_TEST(invoke_to_string(d).size() == expected_size);
        }
    }

    // Day
    for (std::uint8_t day = 0u; day <= 31u; ++day)
    {
        BOOST_TEST_CONTEXT(day)
        {
            datetime d(2021, 1, day);
            BOOST_TEST(invoke_to_string(d).size() == expected_size);
        }
    }

    // Hour
    for (std::uint8_t hour = 0u; hour <= 23u; ++hour)
    {
        BOOST_TEST_CONTEXT(hour)
        {
            datetime d(2021, 1, 3, hour, 10, 15);
            BOOST_TEST(invoke_to_string(d).size() == expected_size);
        }
    }

    // Minute
    for (std::uint8_t minute = 0u; minute <= 59u; ++minute)
    {
        BOOST_TEST_CONTEXT(minute)
        {
            datetime d(2021, 1, 3, 10, minute, 15);
            BOOST_TEST(invoke_to_string(d).size() == expected_size);
        }
    }

    // Second
    for (std::uint8_t second = 0u; second <= 59u; ++second)
    {
        BOOST_TEST_CONTEXT(second)
        {
            datetime d(2021, 1, 3, 10, 43, second);
            BOOST_TEST(invoke_to_string(d).size() == expected_size);
        }
    }

    // Microsecond. Same as for year
    int_generator<std::uint32_t> micro_gen(std::uint32_t(0), std::uint32_t(999999));
    for (int i = 0; i < 50; ++i)
    {
        std::uint32_t micro = micro_gen.generate();
        BOOST_TEST_CONTEXT(micro)
        {
            datetime d(2021, 1, 3, 10, 43, 10, micro);
            BOOST_TEST(invoke_to_string(d).size() == expected_size);
        }
    }
}

// time
static std::string invoke_to_string(::boost::mysql::time t)
{
    std::string res(64, '\0');
    std::size_t sz = detail::time_to_string(t, boost::span<char, 64>(&res[0], 64));
    res.resize(sz);
    return res;
}

BOOST_AUTO_TEST_CASE(time_minmax)
{
    // Double-check C++ min/max work (regression check)
    BOOST_TEST(invoke_to_string((::boost::mysql::time::min)()) == "-2562047788:00:54.775808");
    BOOST_TEST(
        invoke_to_string((::boost::mysql::time::min)() + std::chrono::microseconds(1)) ==
        "-2562047788:00:54.775807"
    );
    BOOST_TEST(invoke_to_string((::boost::mysql::time::max)()) == "2562047788:00:54.775807");
    BOOST_TEST(
        invoke_to_string((::boost::mysql::time::max)() - std::chrono::microseconds(1)) ==
        "2562047788:00:54.775806"
    );
}

BOOST_AUTO_TEST_CASE(time_coverage)
{
    // We use a dot product to cover common cases
    constexpr component_value<int> sign_values[] = {
        {"positive", 1,  "" },
        {"negative", -1, "-"}
    };

    constexpr component_value<int> hours_values[] = {
        {"zero",      0,   "00" },
        {"onedigit",  5,   "05" },
        {"twodigits", 23,  "23" },
        {"max",       838, "838"}
    };

    constexpr component_value<int> mins_secs_values[] = {
        {"zero",      0,  "00"},
        {"onedigit",  5,  "05"},
        {"twodigits", 59, "59"}
    };

    constexpr component_value<int> micros_values[] = {
        {"zero",      0,      "000000"},
        {"onedigit",  5,      "000005"},
        {"twodigits", 50,     "000050"},
        {"max",       999999, "999999"},
    };

    // clang-format off
    for (const auto& sign : sign_values)
    {
    for (const auto& hours : hours_values)
    {
    for (const auto& mins : mins_secs_values)
    {
    for (const auto& secs : mins_secs_values)
    {
    for (const auto& micros : micros_values)
    {
        BOOST_TEST_CONTEXT(
            "sign=" << sign.name << ", hours=" << hours.name <<  ", mins=" << mins.name <<
            ", secs=" << secs.name << "micros=" << micros.name
        )
        {
            // Input value
            auto t = sign.value * test::maket(hours.value, mins.value, secs.value, micros.value);
            if (sign.value == -1 && t == test::maket(0, 0, 0))
            {
                // This case makes no sense, as negative zero is represented as zero
                continue;
            }

            // Expected value
            auto expected = test::stringize(
                sign.repr, hours.repr, ':', mins.repr, ':', secs.repr, '.', micros.repr
            );

            // Call the function
            auto actual = invoke_to_string(t);

            // Check
            BOOST_TEST(actual == expected);

        }
    }
    }
    }
    }
    }
    // clang-format on
}

BOOST_AUTO_TEST_CASE(time_padding)
{
    // Double-check we correctly pad, regardless of the number
    // All times below xx:xx:xx.xxxxxx should have 15 characters
    constexpr std::size_t expected_size = 15;

    // Hour
    for (std::uint8_t hour = 0u; hour <= 99u; ++hour)
    {
        BOOST_TEST_CONTEXT(hour)
        {
            auto t = maket(hour, 11, 20);
            BOOST_TEST(invoke_to_string(t).size() == expected_size);
        }
    }

    // Minute
    for (std::uint8_t minute = 0u; minute <= 59u; ++minute)
    {
        BOOST_TEST_CONTEXT(minute)
        {
            auto t = maket(12, minute, 20);
            BOOST_TEST(invoke_to_string(t).size() == expected_size);
        }
    }

    // Second
    for (std::uint8_t second = 0u; second <= 59u; ++second)
    {
        BOOST_TEST_CONTEXT(second)
        {
            auto t = maket(12, 10, second);
            BOOST_TEST(invoke_to_string(t).size() == expected_size);
        }
    }

    // Microsecond. Iterating over all micros is too costly, so we test some random ones
    int_generator<std::uint32_t> micro_gen(std::uint32_t(0), std::uint32_t(999999));
    for (int i = 0; i < 50; ++i)
    {
        std::uint32_t micro = micro_gen.generate();
        BOOST_TEST_CONTEXT(micro)
        {
            auto t = maket(12, 10, 34, micro);
            BOOST_TEST(invoke_to_string(t).size() == expected_size);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
