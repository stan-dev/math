// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/charconv/issues/166

#include <boost/charconv.hpp>
#include <boost/core/lightweight_test.hpp>
#include <string>

template <typename T>
void test()
{
    constexpr T value = 3746.348756384763L;
    constexpr int precision = 6;

    char buffer[1024];
    const auto result = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), value, boost::charconv::chars_format::fixed, precision);
    *result.ptr = '\0';
    BOOST_TEST(result.ec == std::errc());
    BOOST_TEST_EQ(std::string{buffer}, std::to_string(3746.348756));
}

template <typename T>
void rounding()
{
    constexpr T value = 3746.348759784763L;
    constexpr int precision = 6;

    char buffer[1024];
    const auto result = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), value, boost::charconv::chars_format::fixed, precision);
    *result.ptr = '\0';
    BOOST_TEST(result.ec == std::errc());
    BOOST_TEST_EQ(std::string{buffer}, std::to_string(3746.348760));
}

template <typename T>
void more_rounding()
{
    constexpr T value = 3746.89999999999999999L;
    constexpr int precision = 6;

    char buffer[1024];
    const auto result = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), value, boost::charconv::chars_format::fixed, precision);
    *result.ptr = '\0';
    BOOST_TEST(result.ec == std::errc());
    BOOST_TEST_EQ(std::string{buffer}, std::to_string(3746.900000));
}

template <typename T>
void further_rounding()
{
    const std::array<std::string, 8U> solutions = {{
        "3331.5",
        "3331.52",
        "3331.520",
        "3331.5199",
        "3331.51989",
        "3331.519894",
        "3331.5198945",
        "3331.51989448"
    }};

    constexpr T value = 3331.519894481067353L;

    for (int i = 1; i < 9; ++i)
    {
        char buffer[1024];
        const auto result = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), value, boost::charconv::chars_format::fixed, i);
        *result.ptr = '\0';
        BOOST_TEST(result.ec == std::errc());
        if (!BOOST_TEST_EQ(std::string(buffer), solutions[static_cast<std::size_t>(i - 1)]))
        {
            // LCOV_EXCL_START
            std::cerr << "Precision: " << i
                      << "\nExpected: " << solutions[static_cast<std::size_t>(i - 1)]
                      << "\n     Got: " << buffer << std::endl;
            // LCOV_EXCL_STOP
        }
    }
}

template <typename T>
void full_rounding_test()
{
    constexpr T value = 9999.999999999999999999L;
    constexpr int precision = 6;

    char buffer[1024];
    const auto result = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), value, boost::charconv::chars_format::fixed, precision);
    *result.ptr = '\0';
    BOOST_TEST(result.ec == std::errc());
    BOOST_TEST_EQ(std::string{buffer}, std::to_string(10000.000000));
}

template <typename T>
void test_zeros_path()
{
    constexpr T value = 123456789.000000L;
    constexpr int precision = 6;

    char buffer[1024];
    const auto result = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), value, boost::charconv::chars_format::fixed, precision);
    *result.ptr = '\0';
    BOOST_TEST(result.ec == std::errc());
    BOOST_TEST_EQ(std::string{buffer}, std::to_string(123456789.000000));
}

int main()
{
    #if BOOST_CHARCONV_LDBL_BITS == 80
    test<long double>();
    rounding<long double>();
    more_rounding<long double>();
    further_rounding<long double>();
    full_rounding_test<long double>();
    test_zeros_path<long double>();
    #endif

    #ifdef BOOST_CHARCONV_HAS_QUADMATH
    test<__float128>();
    rounding<__float128>();
    more_rounding<__float128>();
    further_rounding<__float128>();
    full_rounding_test<__float128>();
    test_zeros_path<__float128>();
    #endif

    #if defined(BOOST_CHARCONV_HAS_STDFLOAT128) && defined(BOOST_CHARCONV_HAS_QUADMATH)
    test<std::float128_t>();
    rounding<std::float128_t>();
    more_rounding<std::float128_t>();
    further_rounding<std::float128_t>();
    full_rounding_test<std::float128_t>();
    test_zeros_path<std::float128_t>();
    #endif

    return boost::report_errors();
}
