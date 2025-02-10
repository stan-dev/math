// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/charconv.hpp>
#include <boost/config.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/core/detail/string_view.hpp>
#include <random>
#include <string>
#include <limits>
#include <iomanip>

#if !defined(BOOST_NO_CXX17_HDR_STRING_VIEW)
#  include <string_view>
#endif

static std::mt19937_64 rng(42);
constexpr std::size_t N = 1024;

template <typename T, typename StringViewType = boost::core::string_view>
void test_int()
{
    std::uniform_int_distribution<T> dist((std::numeric_limits<T>::min)(), (std::numeric_limits<T>::max)());

    for (std::size_t i = 0; i < N; ++i)
    {
        const auto value = dist(rng);
        std::string str_value = std::to_string(value);
        StringViewType sv = str_value;

        T v = 0;
        auto r = boost::charconv::from_chars(sv, v);
        BOOST_TEST(r);
        BOOST_TEST_EQ(v, value);
    }
}

template <typename T, typename StringViewType = boost::core::string_view>
void test_float()
{
    std::uniform_real_distribution<T> dist(T(-1e3), T(1e3));

    for (std::size_t i = 0; i < N; ++i)
    {
        const auto value = dist(rng);
        std::string str_value = std::to_string(value);
        StringViewType sv = str_value;

        T v = 0;
        auto r = boost::charconv::from_chars(sv, v);

        T v2 = 0;
        auto r2 = boost::charconv::from_chars(str_value.data(), str_value.data() + str_value.size(), v2);

        if( BOOST_TEST( r.ec == std::errc() ) && BOOST_TEST( r2.ec == std::errc() ) && BOOST_TEST( v == v2 ) )
        {
        }
        else
        {
            std::cerr << std::setprecision(std::numeric_limits<T>::digits10)
                      << "... test failure for value=" << value << "; buffer='" << str_value << "'; errc=" << static_cast<int>(r.ec) << "," << static_cast<int>(r2.ec) << std::endl; // LCOV_EXCL_LINE
        }
    }
}

int main()
{
    // MSVC does not allow (un)signed char in uniform_int_distribution
    // GCC and clang do
    #ifndef _MSC_VER
    test_int<signed char>();
    test_int<unsigned char>();
    #endif

    test_int<short>();
    test_int<unsigned short>();
    test_int<int>();
    test_int<unsigned>();
    test_int<long>();
    test_int<unsigned long>();
    test_int<long long>();
    test_int<unsigned long long>();

    #ifndef _MSC_VER
    test_int<signed char, std::string>();
    test_int<unsigned char, std::string>();
    #endif

    test_int<short, std::string>();
    test_int<unsigned short, std::string>();
    test_int<int, std::string>();
    test_int<unsigned, std::string>();
    test_int<long, std::string>();
    test_int<unsigned long, std::string>();
    test_int<long long, std::string>();
    test_int<unsigned long long, std::string>();

    #if !defined(BOOST_NO_CXX17_HDR_STRING_VIEW)

    #ifndef _MSC_VER
    test_int<signed char, std::string_view>();
    test_int<unsigned char, std::string_view>();
    #endif

    test_int<short, std::string_view>();
    test_int<unsigned short, std::string_view>();
    test_int<int, std::string_view>();
    test_int<unsigned, std::string_view>();
    test_int<long, std::string_view>();
    test_int<unsigned long, std::string_view>();
    test_int<long long, std::string_view>();
    test_int<unsigned long long, std::string_view>();

    #endif

    test_float<float>();
    test_float<double>();

    test_float<float, std::string>();
    test_float<double, std::string>();

    #ifndef BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE
    test_float<long double>();
    test_float<long double, std::string>();
    #endif

    #if !defined(BOOST_NO_CXX17_HDR_STRING_VIEW)

    test_float<float, std::string_view>();
    test_float<double, std::string_view>();

    #ifndef BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE
    test_float<long double, std::string_view>();
    #endif

    #endif

    #ifdef BOOST_CHARCONV_HAS_FLOAT32

    test_float<std::float32_t>();
    test_float<std::float32_t, std::string>();
    test_float<std::float32_t, std::string_view>();

    #endif

    #ifdef BOOST_CHARCONV_HAS_FLOAT64

    test_float<std::float64_t>();
    test_float<std::float64_t, std::string>();
    test_float<std::float64_t, std::string_view>();

    #endif

    return boost::report_errors();
}
