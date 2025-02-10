// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/charconv/detail/parser.hpp>
#include <boost/charconv/chars_format.hpp>
#include <boost/core/lightweight_test.hpp>
#include <system_error>
#include <type_traits>
#include <cstdint>
#include <cstring>
#include <cerrno>

template <typename T>
void test_integer()
{
    std::uint64_t significand {};
    std::int64_t  exponent {};
    bool sign {};

    const char* val1 = "12";
    auto r1 = boost::charconv::detail::parser(val1, val1 + std::strlen(val1), sign, significand, exponent);
    BOOST_TEST(r1.ec == std::errc());
    BOOST_TEST_EQ(sign, false);
    BOOST_TEST_EQ(significand, UINT64_C(12));
    BOOST_TEST_EQ(exponent, 0);

    significand = 0;
    exponent = 0;
    sign = false;

    const char* val2 = "123456789";
    auto r2 = boost::charconv::detail::parser(val2, val2 + std::strlen(val2), sign, significand, exponent);
    BOOST_TEST(r2.ec == std::errc());
    BOOST_TEST_EQ(sign, false);
    BOOST_TEST_EQ(exponent, 0);
    BOOST_TEST_EQ(significand, UINT64_C(123456789));

    auto r3 = boost::charconv::detail::parser(val2, val2 + std::strlen(val2), sign, significand, exponent, boost::charconv::chars_format::scientific);
    BOOST_TEST(r3.ec == std::errc::invalid_argument);
}

template <typename T>
void test_scientifc()
{
    std::uint64_t significand {};
    std::int64_t  exponent {};
    bool sign {};

    const char* val1 = "-1e1";
    auto r1 = boost::charconv::detail::parser(val1, val1 + std::strlen(val1), sign, significand, exponent);
    BOOST_TEST(r1.ec == std::errc());
    BOOST_TEST_EQ(sign, true);
    BOOST_TEST_EQ(significand, UINT64_C(1));
    BOOST_TEST_EQ(exponent, 1);

    significand = 0;
    exponent = 0;
    sign = false;

    const char* val2 = "123456789e10";
    auto r2 = boost::charconv::detail::parser(val2, val2 + std::strlen(val2), sign, significand, exponent);
    BOOST_TEST(r2.ec == std::errc());
    BOOST_TEST_EQ(sign, false);
    BOOST_TEST_EQ(exponent, 10);
    BOOST_TEST_EQ(significand, UINT64_C(123456789));

    significand = 0;
    exponent = 0;
    sign = false;

    const char* val3 = "1.23456789e+10";
    auto r3 = boost::charconv::detail::parser(val3, val3 + std::strlen(val3), sign, significand, exponent);
    BOOST_TEST(r3.ec == std::errc());
    BOOST_TEST_EQ(sign, false);
    BOOST_TEST_EQ(exponent, 2);
    BOOST_TEST_EQ(significand, UINT64_C(123456789));

    const char* val4 = "1.23456789e-10";
    auto r4 = boost::charconv::detail::parser(val4, val4 + std::strlen(val4), sign, significand, exponent);
    BOOST_TEST(r4.ec == std::errc());
    BOOST_TEST_EQ(sign, false);
    BOOST_TEST_EQ(exponent, -18);
    BOOST_TEST_EQ(significand, UINT64_C(123456789));

    auto r5 = boost::charconv::detail::parser(val4, val4 + std::strlen(val4), sign, significand, exponent, boost::charconv::chars_format::fixed);
    BOOST_TEST(r5.ec == std::errc::invalid_argument);

    const char* val6 = "987654321e10";
    auto r6 = boost::charconv::detail::parser(val6, val6 + std::strlen(val6), sign, significand, exponent);
    BOOST_TEST(r6.ec == std::errc());
    BOOST_TEST_EQ(sign, false);
    BOOST_TEST_EQ(exponent, 10);
    BOOST_TEST_EQ(significand, UINT64_C(987654321));

    const char* val7 = "1.23456789E+10";
    auto r7 = boost::charconv::detail::parser(val7, val7 + std::strlen(val7), sign, significand, exponent);
    BOOST_TEST(r7.ec == std::errc());
    BOOST_TEST_EQ(sign, false);
    BOOST_TEST_EQ(exponent, 2);
    BOOST_TEST_EQ(significand, UINT64_C(123456789));

}

template <typename T>
void test_hex_integer()
{
    std::uint64_t significand {};
    std::int64_t  exponent {};
    bool sign {};

    const char* val1 = "2a";
    auto r1 = boost::charconv::detail::parser(val1, val1 + std::strlen(val1), sign, significand, exponent, boost::charconv::chars_format::hex);
    BOOST_TEST(r1.ec == std::errc());
    BOOST_TEST_EQ(sign, false);
    BOOST_TEST_EQ(significand, UINT64_C(42));
    BOOST_TEST_EQ(exponent, 0);

    significand = 0;
    exponent = 0;
    sign = false;

    const char* val2 = "-1a3b5c7d9";
    auto r2 = boost::charconv::detail::parser(val2, val2 + std::strlen(val2), sign, significand, exponent, boost::charconv::chars_format::hex);
    BOOST_TEST(r2.ec == std::errc());
    BOOST_TEST_EQ(sign, true);
    BOOST_TEST_EQ(exponent, 0);
    BOOST_TEST_EQ(significand, UINT64_C(7041566681));

    auto r3 = boost::charconv::detail::parser(val2, val2 + std::strlen(val2), sign, significand, exponent, boost::charconv::chars_format::scientific);
    BOOST_TEST(r3.ec == std::errc::invalid_argument);
}

template <typename T>
void test_hex_scientific()
{
    std::uint64_t significand {};
    std::int64_t  exponent {};
    bool sign {};

    const char* val1 = "2ap+5";
    auto r1 = boost::charconv::detail::parser(val1, val1 + std::strlen(val1), sign, significand, exponent, boost::charconv::chars_format::hex);
    BOOST_TEST(r1.ec == std::errc());
    BOOST_TEST_EQ(sign, false);
    BOOST_TEST_EQ(significand, UINT64_C(42));
    BOOST_TEST_EQ(exponent, 5);

    significand = 0;
    exponent = 0;
    sign = false;

    const char* val2 = "-1.3a2bp-10";
    auto r2 = boost::charconv::detail::parser(val2, val2 + std::strlen(val2), sign, significand, exponent, boost::charconv::chars_format::hex);
    BOOST_TEST(r2.ec == std::errc());
    BOOST_TEST_EQ(sign, true);
    BOOST_TEST_EQ(exponent, -14);
    BOOST_TEST_EQ(significand, UINT64_C(80427));

    auto r3 = boost::charconv::detail::parser(val2, val2 + std::strlen(val2), sign, significand, exponent, boost::charconv::chars_format::scientific);
    BOOST_TEST(r3.ec == std::errc::invalid_argument);

    const char* val4 = "-1.3A2BP-10";
    auto r4 = boost::charconv::detail::parser(val4, val4 + std::strlen(val4), sign, significand, exponent, boost::charconv::chars_format::hex);
    BOOST_TEST(r4.ec == std::errc());
    BOOST_TEST_EQ(sign, true);
    BOOST_TEST_EQ(exponent, -14);
    BOOST_TEST_EQ(significand, UINT64_C(80427));
}

int main()
{
    test_integer<float>();
    test_integer<double>();
    test_integer<long double>();

    test_scientifc<float>();
    test_scientifc<double>();
    test_scientifc<long double>();

    test_hex_integer<float>();
    test_hex_integer<double>();
    test_hex_integer<long double>();

    test_hex_scientific<float>();
    test_hex_scientific<double>();
    test_hex_scientific<long double>();

    return boost::report_errors();
}
