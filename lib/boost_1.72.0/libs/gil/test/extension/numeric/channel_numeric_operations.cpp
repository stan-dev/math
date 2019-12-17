//
// Copyright 2019 Mateusz Loskot <mateusz at loskot dot net>
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//
#include <boost/gil.hpp>
#include <boost/gil/extension/numeric/channel_numeric_operations.hpp>

#include <tuple>
#include <type_traits>

#define BOOST_TEST_MODULE test_ext_numeric_pixel_numeric_operations
#include "unit_test.hpp"
#include "unit_test_utility.hpp"
#include "core/channel/test_fixture.hpp"

namespace gil = boost::gil;
namespace fixture = boost::gil::test::fixture;

BOOST_AUTO_TEST_SUITE(channel_plus_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(plus_integer_same_types, channel_t, fixture::channel_integer_types)
{
    gil::channel_plus_t<channel_t, channel_t, channel_t> f;
    BOOST_TEST(f(0, 0) == channel_t(0));
    BOOST_TEST(f(100, 27) == channel_t(127));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(plus_integer_mixed_types, channel_t, fixture::channel_integer_types)
{
    {
        using channel1_t = channel_t;
        using channel2_t = std::uint8_t; // duplicates only one of fixture::channel_integer_types
        gil::channel_plus_t<channel1_t, channel2_t, channel1_t> f;
        BOOST_TEST(f(0, 0) == channel1_t(0));
        BOOST_TEST(f(100, 27) == channel_t(127));
    }
    {
        using channel1_t = std::uint8_t; // duplicates only one of fixture::channel_integer_types
        using channel2_t = channel_t;
        gil::channel_plus_t<channel1_t, channel2_t, channel2_t> f;
        BOOST_TEST(f(0, 0) == channel2_t(0));
        BOOST_TEST(f(100, 27) == channel_t(127));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(plus_integer_signed_types_with_overflow, channel_t, fixture::channel_integer_signed_types)
{
    // Signed integer overflow is UB, so just check addition does not yield mathematically
    // expected value but is constrained by the range of representable values for given type.

    auto const max_value = gil::channel_traits<channel_t>::max_value();
    gil::channel_plus_t<channel_t, channel_t, channel_t> f;
    BOOST_TEST(f(max_value, 1) != std::int64_t(max_value) + 1);
    BOOST_TEST(f(max_value, max_value) != std::int64_t(max_value) + max_value);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(plus_integer_unsigned_types_with_wraparound, channel_t, fixture::channel_integer_unsigned_types)
{
    // The C Standard, 6.2.5, paragraph 9 [ISO/IEC 9899:2011], states:
    // A computation involving unsigned operands can never overflow, because a result that
    // cannot be represented by the resulting unsigned integer type is reduced modulo the number
    // that is one greater than the largest value that can be represented by the resulting type.

    auto const max_value = gil::channel_traits<channel_t>::max_value();
    auto const min_value = gil::channel_traits<channel_t>::min_value();
    gil::channel_plus_t<channel_t, channel_t, channel_t> f;
    BOOST_TEST(f(max_value, 1) == min_value);
    BOOST_TEST(f(max_value, max_value) == max_value - 1);
}

BOOST_AUTO_TEST_SUITE_END() // channel_plus_t

BOOST_AUTO_TEST_SUITE(channel_minus_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(minus_integer_same_types, channel_t, fixture::channel_integer_types)
{
    gil::channel_minus_t<channel_t, channel_t, channel_t> f;
    BOOST_TEST(f(0, 0) == channel_t(0));
    BOOST_TEST(f(100, 27) == channel_t(73));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(minus_integer_mixed_types, channel_t, fixture::channel_integer_types)
{
    {
        using channel1_t = channel_t;
        using channel2_t = std::uint8_t; // duplicates only one of fixture::channel_integer_types
        gil::channel_minus_t<channel1_t, channel2_t, channel1_t> f;
        BOOST_TEST(f(0, 0) == channel1_t(0));
        BOOST_TEST(f(100, 27) == channel_t(73));
    }
    {
        using channel1_t = std::uint8_t; // duplicates only one of fixture::channel_integer_types
        using channel2_t = channel_t;
        gil::channel_minus_t<channel1_t, channel2_t, channel2_t> f;
        BOOST_TEST(f(0, 0) == channel2_t(0));
        BOOST_TEST(f(100, 27) == channel_t(73));
    }
}

BOOST_AUTO_TEST_SUITE_END() // channel_minus_t

BOOST_AUTO_TEST_SUITE(channel_multiplies_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(multiplies_integer_same_types, channel_t, fixture::channel_integer_types)
{
    gil::channel_multiplies_t<channel_t, channel_t, channel_t> f;
    BOOST_TEST(f(0, 0) == channel_t(0));
    BOOST_TEST(f(1, 1) == channel_t(1));
    BOOST_TEST(f(4, 2) == channel_t(8));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multiplies_integer_mixed_types, channel_t, fixture::channel_integer_types)
{
    {
        using channel1_t = channel_t;
        using channel2_t = std::uint8_t; // duplicates only one of fixture::channel_integer_types
        gil::channel_multiplies_t<channel1_t, channel2_t, channel1_t> f;
        BOOST_TEST(f(0, 0) == channel1_t(0));
        BOOST_TEST(f(4, 2) == channel_t(8));
    }
    {
        using channel1_t = std::uint8_t; // duplicates only one of fixture::channel_integer_types
        using channel2_t = channel_t;
        gil::channel_multiplies_t<channel1_t, channel2_t, channel2_t> f;
        BOOST_TEST(f(0, 0) == channel2_t(0));
        BOOST_TEST(f(4, 2) == channel_t(8));
    }
}

BOOST_AUTO_TEST_SUITE_END() // channel_multiplies_t

BOOST_AUTO_TEST_SUITE(channel_divides_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(divides_integer_same_types, channel_t, fixture::channel_integer_types)
{
    gil::channel_divides_t<channel_t, channel_t, channel_t> f;
    BOOST_TEST(f(0, 1) == channel_t(0));
    BOOST_TEST(f(1, 1) == channel_t(1));
    BOOST_TEST(f(4, 2) == channel_t(2));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(divides_integer_mixed_types, channel_t, fixture::channel_integer_types)
{
    {
        using channel1_t = channel_t;
        using channel2_t = std::uint8_t; // duplicates only one of fixture::channel_integer_types
        gil::channel_divides_t<channel1_t, channel2_t, channel1_t> f;
        BOOST_TEST(f(0, 1) == channel1_t(0));
        BOOST_TEST(f(4, 2) == channel_t(2));
    }
    {
        using channel1_t = std::uint8_t; // duplicates only one of fixture::channel_integer_types
        using channel2_t = channel_t;
        gil::channel_divides_t<channel1_t, channel2_t, channel2_t> f;
        BOOST_TEST(f(0, 1) == channel2_t(0));
        BOOST_TEST(f(4, 2) == channel_t(2));
    }
}

BOOST_AUTO_TEST_SUITE_END() // channel_divides_t

BOOST_AUTO_TEST_SUITE(channel_plus_scalar_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(plus_scalar_integer_same_types, channel_t, fixture::channel_integer_types)
{
    gil::channel_plus_scalar_t<channel_t, int, channel_t> f;
    BOOST_TEST(f(0, 0) == channel_t(0));
    BOOST_TEST(f(100, 27) == channel_t(127));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(plus_scalar_integer_mixed_types, channel_t, fixture::channel_integer_types)
{
    using channel_result_t = std::uint8_t;
    gil::channel_plus_scalar_t<channel_t, int, channel_result_t> f;
    BOOST_TEST(f(0, 0) == channel_result_t(0));
    BOOST_TEST(f(100, 27) == channel_result_t(127));
}

BOOST_AUTO_TEST_SUITE_END() // channel_plus_scalar_t

BOOST_AUTO_TEST_SUITE(channel_minus_scalar_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(minus_scalar_integer_same_types, channel_t, fixture::channel_integer_types)
{
    gil::channel_minus_scalar_t<channel_t, int, channel_t> f;
    BOOST_TEST(f(0, 0) == channel_t(0));
    BOOST_TEST(f(100, 27) == channel_t(73));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(minus_scalar_integer_mixed_types, channel_t, fixture::channel_integer_types)
{
    using channel_result_t = std::uint8_t;
    gil::channel_minus_scalar_t<channel_t, int, std::uint8_t> f;
    BOOST_TEST(f(0, 0) == channel_result_t(0));
    BOOST_TEST(f(100, 27) == channel_result_t(73));
}

BOOST_AUTO_TEST_SUITE_END() // channel_minus_scalar_t

BOOST_AUTO_TEST_SUITE(channel_multiplies_scalar_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(multiplies_scalar_integer_same_types, channel_t, fixture::channel_integer_types)
{
    gil::channel_multiplies_scalar_t<channel_t, channel_t, channel_t> f;
    BOOST_TEST(f(0, 0) == channel_t(0));
    BOOST_TEST(f(1, 1) == channel_t(1));
    BOOST_TEST(f(4, 2) == channel_t(8));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multiplies_scalar_integer_mixed_types, channel_t, fixture::channel_integer_types)
{
    using channel_result_t = std::uint8_t;
    gil::channel_multiplies_scalar_t<channel_t, int, channel_result_t> f;
    BOOST_TEST(f(0, 0) == channel_result_t(0));
    BOOST_TEST(f(4, 2) == channel_result_t(8));
}

BOOST_AUTO_TEST_SUITE_END() // channel_multiplies_scalar_t

BOOST_AUTO_TEST_SUITE(channel_divides_scalar_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(divides_scalar_integer_same_types, channel_t, fixture::channel_integer_types)
{
    gil::channel_divides_scalar_t<channel_t, channel_t, channel_t> f;
    BOOST_TEST(f(0, 1) == channel_t(0));
    BOOST_TEST(f(1, 1) == channel_t(1));
    BOOST_TEST(f(4, 2) == channel_t(2));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(divides_scalar_integer_mixed_types, channel_t, fixture::channel_integer_types)
{
    using channel_result_t = std::uint8_t; // duplicates only one of fixture::channel_integer_types
    gil::channel_divides_scalar_t<channel_t, int, channel_result_t> f;
    BOOST_TEST(f(0, 1) == channel_t(0));
    BOOST_TEST(f(4, 2) == channel_t(2));
}

BOOST_AUTO_TEST_SUITE_END() // channel_divides_scalar_t

BOOST_AUTO_TEST_SUITE(channel_halves_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(halves_integer_same_types, channel_t, fixture::channel_integer_types)
{
    gil::channel_halves_t<channel_t> f;
    {
        channel_t c(0);
        f(c);
        BOOST_TEST(c == channel_t(0));
    }
    {
        channel_t c(2);
        f(c);
        BOOST_TEST(c == channel_t(1));
    }
    {
        channel_t c(4);
        f(c);
        BOOST_TEST(c == channel_t(2));
    }
}

BOOST_AUTO_TEST_SUITE_END() // channel_halves_t

BOOST_AUTO_TEST_SUITE(channel_zeros_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(zeros_integer_same_types, channel_t, fixture::channel_integer_types)
{
    gil::channel_zeros_t<channel_t> f;
    {
        channel_t c(0);
        f(c);
        BOOST_TEST(c == channel_t(0));
    }
    {
        channel_t c(2);
        f(c);
        BOOST_TEST(c == channel_t(0));
    }
    {
        channel_t c(4);
        f(c);
        BOOST_TEST(c == channel_t(0));
    }
}

BOOST_AUTO_TEST_SUITE_END() // channel_zeros_t

BOOST_AUTO_TEST_SUITE(channel_assigns_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(assigns_integer_same_types, channel_t, fixture::channel_integer_types)
{
    gil::channel_assigns_t<channel_t, channel_t> f;
    {
        channel_t c1(10);
        channel_t c2(20);
        f(c1, c2);
        BOOST_TEST(c2 == c1);
    }

}

BOOST_AUTO_TEST_SUITE_END() // channel_assigns_t
