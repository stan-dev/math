// Copyright 2018 Ulf Adams
// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/charconv.hpp>
#include <boost/charconv/detail/fallback_routines.hpp>
#include <boost/core/lightweight_test.hpp>
#include <system_error>
#include <type_traits>
#include <limits>
#include <cstring>
#include <cstdint>
#include <cerrno>
#include <utility>
#include <string>
#include <random>
#include <iomanip>
#include <sstream>
#include <boost/core/detail/splitmix64.hpp>

// These numbers diverge from what the formatting is using printf
// See: https://godbolt.org/z/zd34KcWMW
template <typename T>
void printf_divergence()
{
    char buffer1[256] {};
    T v1 = 3.4;
    auto r1 = boost::charconv::to_chars(buffer1, buffer1 + sizeof(buffer1), v1);
    BOOST_TEST(r1.ec == std::errc());
    BOOST_TEST_CSTR_EQ(buffer1, "3.4");

    char buffer2[256] {};
    T v2 = 3000.40;
    auto r2 = boost::charconv::to_chars(buffer2, buffer2 + sizeof(buffer2), v2);
    BOOST_TEST(r2.ec == std::errc());
    BOOST_TEST_CSTR_EQ(buffer2, "3000.4");

    char buffer3[256] {};
    T v3 = -3000000300000000.5;
    auto r3 = boost::charconv::to_chars(buffer3, buffer3 + sizeof(buffer3), v3);
    BOOST_TEST(r3.ec == std::errc());
    BOOST_TEST_CSTR_EQ(buffer3, "-3000000300000000.5");
}

template <typename T>
void integer_general_format()
{
    char buffer1[256] {};
    T v1 = 1217.2772861138403;
    auto r1 = boost::charconv::to_chars(buffer1, buffer1 + sizeof(buffer1), v1);
    BOOST_TEST(r1.ec == std::errc());
    BOOST_TEST_CSTR_EQ(buffer1, "1217.2772861138403");
    T return_v1;
    auto r1_return = boost::charconv::from_chars(buffer1, buffer1 + strlen(buffer1), return_v1);
    BOOST_TEST(r1_return.ec == std::errc());
    BOOST_TEST_EQ(return_v1, v1);
}

template <typename T>
void non_finite_values(boost::charconv::chars_format fmt = boost::charconv::chars_format::general, int precision = -1)
{
    char buffer1[256] {};
    T v1 = std::numeric_limits<T>::infinity();
    auto r1 = boost::charconv::to_chars(buffer1, buffer1 + sizeof(buffer1), v1, fmt, precision);
    BOOST_TEST(r1.ec == std::errc());
    BOOST_TEST_CSTR_EQ(buffer1, "inf");

    char buffer2[256] {};
    T v2 = -std::numeric_limits<T>::infinity();
    auto r2 = boost::charconv::to_chars(buffer2, buffer2 + sizeof(buffer2), v2, fmt, precision);
    BOOST_TEST(r2.ec == std::errc());
    BOOST_TEST_CSTR_EQ(buffer2, "-inf");

    char buffer3[256] {};
    T v3 = std::numeric_limits<T>::quiet_NaN();
    auto r3 = boost::charconv::to_chars(buffer3, buffer3 + sizeof(buffer3), v3, fmt, precision);
    BOOST_TEST(r3.ec == std::errc());
    BOOST_TEST_CSTR_EQ(buffer3, "nan");

    char buffer4[256] {};
    T v4 = -std::numeric_limits<T>::quiet_NaN();
    auto r4 = boost::charconv::to_chars(buffer4, buffer4 + sizeof(buffer4), v4, fmt, precision);
    BOOST_TEST(r4.ec == std::errc());
    BOOST_TEST_CSTR_EQ(buffer4, "-nan(ind)");

    char buffer5[256] {};
    T v5 = std::numeric_limits<T>::signaling_NaN();
    auto r5 = boost::charconv::to_chars(buffer5, buffer5 + sizeof(buffer5), v5, fmt, precision);
    BOOST_TEST(r5.ec == std::errc());
    BOOST_TEST_CSTR_EQ(buffer5, "nan(snan)");

    char buffer6[256] {};
    T v6 = -std::numeric_limits<T>::signaling_NaN();
    auto r6 = boost::charconv::to_chars(buffer6, buffer6 + sizeof(buffer6), v6, fmt, precision);
    BOOST_TEST(r6.ec == std::errc());
    BOOST_TEST_CSTR_EQ(buffer6, "-nan(snan)");
}

template <typename T>
void fixed_values()
{
    char buffer1[256] {};
    T v1 = 61851632;
    auto r1 = boost::charconv::to_chars(buffer1, buffer1 + sizeof(buffer1), v1);
    BOOST_TEST(r1.ec == std::errc());
    BOOST_TEST_CSTR_EQ(buffer1, "61851632");
}

template <typename T>
void failing_ci_values()
{
    char buffer1[256] {};
    T v1 = -1.08260383390082946e+307;
    auto r1 = boost::charconv::to_chars(buffer1, buffer1 + sizeof(buffer1), v1, boost::charconv::chars_format::hex);
    BOOST_TEST(r1.ec == std::errc());
    BOOST_TEST_CSTR_EQ(buffer1, "-1.ed5658af91a0fp+1019");

    char buffer2[256] {};
    T v2 = -9.52743282403084637e+306;
    auto r2 = boost::charconv::to_chars(buffer2, buffer2 + sizeof(buffer2), v2, boost::charconv::chars_format::hex);
    BOOST_TEST(r2.ec == std::errc());
    BOOST_TEST_CSTR_EQ(buffer2, "-1.b22914956c56fp+1019");
}

template <typename T>
void spot_check(T v, const std::string& str, boost::charconv::chars_format fmt)
{
    char buffer[256] {};
    const auto r = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), v, fmt);
    BOOST_TEST(r.ec == std::errc());
    BOOST_TEST_CSTR_EQ(buffer, str.c_str());
}

template <typename T>
void spot_check(T v, const std::string& str, boost::charconv::chars_format fmt, int precision)
{
    char buffer[256] {};
    const auto r = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), v, fmt, precision);
    BOOST_TEST(r.ec == std::errc());
    BOOST_TEST_CSTR_EQ(buffer, str.c_str());
}

constexpr const char* fmt_general(double)
{
    return "%g";
}

constexpr const char* fmt_general(long double)
{
    return "%Lg";
}

constexpr const char* fmt_sci(double)
{
    return "%e";
}

constexpr const char* fmt_sci(long double)
{
    return "%Le";
}

constexpr const char* fmt_fixed(double)
{
    return "%.0f";
}

constexpr const char* fmt_fixed(long double)
{
    return "%.0Lf";
}

template <typename T>
void test_printf_fallback(T v, const std::string&, boost::charconv::chars_format fmt = boost::charconv::chars_format::general, int precision = -1)
{
    char buffer[256] {};
    const auto r = boost::charconv::detail::to_chars_printf_impl(buffer, buffer + sizeof(buffer), v, fmt, precision);
    BOOST_TEST(r.ec == std::errc());

    char printf_buffer[256] {};
    if (fmt == boost::charconv::chars_format::general)
    {
        std::snprintf(printf_buffer, sizeof(printf_buffer), fmt_general(v), v);
    }
    else if (fmt == boost::charconv::chars_format::scientific)
    {
        std::snprintf(printf_buffer, sizeof(printf_buffer), fmt_sci(v), v);
    }
    else if (fmt == boost::charconv::chars_format::fixed)
    {
        std::snprintf(printf_buffer, sizeof(printf_buffer), fmt_fixed(v), v);
    }

    BOOST_TEST_CSTR_EQ(buffer, printf_buffer);
}

std::string format(int prec)
{
    std::string format = "%." + std::to_string(prec) + "g";
    return format;
}

template <typename T>
void test_to_chars_hex_16(T v)
{
    if (!std::isnormal(v))
    {
        return;
    }

    char buffer[256];
    // Stringstream will write with leading 0x whereas we will not
    auto ptr = buffer;
    *ptr++ = '0';
    *ptr++ = 'x';

    auto r = boost::charconv::to_chars(ptr, ptr + sizeof(buffer) - 2, v, boost::charconv::chars_format::hex);
    BOOST_TEST(r);
    *r.ptr = '\0';

    std::stringstream o_val;
    o_val << std::hexfloat << v;

    BOOST_TEST_CSTR_EQ(buffer, o_val.str().c_str());
}

int main()
{
    printf_divergence<double>();
    integer_general_format<double>();
    non_finite_values<double>();
    non_finite_values<double>(boost::charconv::chars_format::general, 2);
    non_finite_values<double>(boost::charconv::chars_format::scientific);
    non_finite_values<double>(boost::charconv::chars_format::scientific, 2);
    non_finite_values<double>(boost::charconv::chars_format::hex);
    non_finite_values<double>(boost::charconv::chars_format::hex, 2);
    non_finite_values<double>(boost::charconv::chars_format::fixed);
    non_finite_values<double>(boost::charconv::chars_format::fixed, 2);

    // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=57484
    #if !(defined(__GNUC__) && __GNUC__ == 4 && __GNUC_MINOR__ < 9 && defined(__i686__)) && !defined(BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE)
    non_finite_values<long double>();
    #endif

    fixed_values<float>();
    fixed_values<double>();

    failing_ci_values<double>();

    // Values from ryu tests
    spot_check(1.0, "1", boost::charconv::chars_format::general);
    spot_check(1.2, "1.2", boost::charconv::chars_format::general);
    spot_check(1.23, "1.23", boost::charconv::chars_format::general);
    spot_check(1.234, "1.234", boost::charconv::chars_format::general);
    spot_check(1.2345, "1.2345", boost::charconv::chars_format::general);
    spot_check(1.23456, "1.23456", boost::charconv::chars_format::general);
    spot_check(1.234567, "1.234567", boost::charconv::chars_format::general);
    spot_check(1.2345678, "1.2345678", boost::charconv::chars_format::general);
    spot_check(1.23456789, "1.23456789", boost::charconv::chars_format::general);
    spot_check(1.234567890, "1.23456789", boost::charconv::chars_format::general);
    spot_check(1.2345678901, "1.2345678901", boost::charconv::chars_format::general);
    spot_check(1.23456789012, "1.23456789012", boost::charconv::chars_format::general);
    spot_check(1.234567890123, "1.234567890123", boost::charconv::chars_format::general);
    spot_check(1.2345678901234, "1.2345678901234", boost::charconv::chars_format::general);
    spot_check(1.23456789012345, "1.23456789012345", boost::charconv::chars_format::general);
    spot_check(1.234567890123456, "1.234567890123456", boost::charconv::chars_format::general);

    spot_check(1.0, "1e+00", boost::charconv::chars_format::scientific);
    spot_check(1.2, "1.2e+00", boost::charconv::chars_format::scientific);
    spot_check(1.23, "1.23e+00", boost::charconv::chars_format::scientific);
    spot_check(1.234, "1.234e+00", boost::charconv::chars_format::scientific);
    spot_check(1.2345, "1.2345e+00", boost::charconv::chars_format::scientific);
    spot_check(1.23456, "1.23456e+00", boost::charconv::chars_format::scientific);
    spot_check(1.234567, "1.234567e+00", boost::charconv::chars_format::scientific);
    spot_check(1.2345678, "1.2345678e+00", boost::charconv::chars_format::scientific);
    spot_check(1.23456789, "1.23456789e+00", boost::charconv::chars_format::scientific);
    spot_check(1.234567890, "1.23456789e+00", boost::charconv::chars_format::scientific);
    spot_check(1.2345678901, "1.2345678901e+00", boost::charconv::chars_format::scientific);
    spot_check(1.23456789012, "1.23456789012e+00", boost::charconv::chars_format::scientific);
    spot_check(1.234567890123, "1.234567890123e+00", boost::charconv::chars_format::scientific);
    spot_check(1.2345678901234, "1.2345678901234e+00", boost::charconv::chars_format::scientific);
    spot_check(1.23456789012345, "1.23456789012345e+00", boost::charconv::chars_format::scientific);
    spot_check(1.234567890123456, "1.234567890123456e+00", boost::charconv::chars_format::scientific);

    spot_check(1.0, "1e+00", boost::charconv::chars_format::scientific);
    spot_check(12.0, "1.2e+01", boost::charconv::chars_format::scientific);
    spot_check(123.0, "1.23e+02", boost::charconv::chars_format::scientific);
    spot_check(1234.0, "1.234e+03", boost::charconv::chars_format::scientific);
    spot_check(12345.0, "1.2345e+04", boost::charconv::chars_format::scientific);
    spot_check(123456.0, "1.23456e+05", boost::charconv::chars_format::scientific);
    spot_check(1234567.0, "1.234567e+06", boost::charconv::chars_format::scientific);
    spot_check(12345678.0, "1.2345678e+07", boost::charconv::chars_format::scientific);
    spot_check(123456789.0, "1.23456789e+08", boost::charconv::chars_format::scientific);
    spot_check(1234567890.0, "1.23456789e+09", boost::charconv::chars_format::scientific);
    spot_check(12345678901.0, "1.2345678901e+10", boost::charconv::chars_format::scientific);
    spot_check(123456789012.0, "1.23456789012e+11", boost::charconv::chars_format::scientific);
    spot_check(1234567890123.0, "1.234567890123e+12", boost::charconv::chars_format::scientific);
    spot_check(12345678901234.0, "1.2345678901234e+13", boost::charconv::chars_format::scientific);
    spot_check(123456789012345.0, "1.23456789012345e+14", boost::charconv::chars_format::scientific);
    spot_check(1234567890123456.0, "1.234567890123456e+15", boost::charconv::chars_format::scientific);

    test_printf_fallback(1.0, "1");
    test_printf_fallback(1.2, "1.2");
    test_printf_fallback(1.23, "1.23");
    test_printf_fallback(1.234, "1.234");
    test_printf_fallback(1.2345, "1.2345");
    test_printf_fallback(1.23456, "1.23456");
    test_printf_fallback(1.234567, "1.234567");
    test_printf_fallback(1.2345678, "1.2345678");
    test_printf_fallback(1.23456789, "1.23456789");
    test_printf_fallback(1.234567890, "1.23456789");
    test_printf_fallback(1.2345678901, "1.2345678901");
    test_printf_fallback(1.23456789012, "1.23456789012");
    test_printf_fallback(1.234567890123, "1.234567890123");
    test_printf_fallback(1.2345678901234, "1.2345678901234");
    test_printf_fallback(1.23456789012345, "1.23456789012345");
    test_printf_fallback(1.234567890123456, "1.234567890123456");

    test_printf_fallback(1.0, "1e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.2, "1.2e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.23, "1.23e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.234, "1.234e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.2345, "1.2345e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.23456, "1.23456e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.234567, "1.234567e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.2345678, "1.2345678e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.23456789, "1.23456789e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.234567890, "1.23456789e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.2345678901, "1.2345678901e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.23456789012, "1.23456789012e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.234567890123, "1.234567890123e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.2345678901234, "1.2345678901234e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.23456789012345, "1.23456789012345e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.234567890123456, "1.234567890123456e+00", boost::charconv::chars_format::scientific);

    test_printf_fallback(1.0, "1e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.2, "1.2e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.23, "1.23e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.234, "1.234e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.2345, "1.2345e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.23456, "1.23456e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.234567, "1.234567e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.2345678, "1.2345678e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.23456789, "1.23456789e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.234567890, "1.23456789e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.2345678901, "1.2345678901e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.23456789012, "1.23456789012e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.234567890123, "1.234567890123e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.2345678901234, "1.2345678901234e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.23456789012345, "1.23456789012345e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.234567890123456, "1.234567890123456e+00", boost::charconv::chars_format::fixed);

    test_printf_fallback(1.0L, "1");
    test_printf_fallback(1.2L, "1.2");
    test_printf_fallback(1.23L, "1.23");
    test_printf_fallback(1.234L, "1.234");
    test_printf_fallback(1.2345L, "1.2345");
    test_printf_fallback(1.23456L, "1.23456");
    test_printf_fallback(1.234567L, "1.234567");
    test_printf_fallback(1.2345678L, "1.2345678");
    test_printf_fallback(1.23456789L, "1.23456789");
    test_printf_fallback(1.234567890L, "1.23456789");
    test_printf_fallback(1.2345678901L, "1.2345678901");
    test_printf_fallback(1.23456789012L, "1.23456789012");
    test_printf_fallback(1.234567890123L, "1.234567890123");
    test_printf_fallback(1.2345678901234L, "1.2345678901234");
    test_printf_fallback(1.23456789012345L, "1.23456789012345");
    test_printf_fallback(1.234567890123456L, "1.234567890123456");

    test_printf_fallback(1.0L, "1e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.2L, "1.2e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.23L, "1.23e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.234L, "1.234e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.2345L, "1.2345e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.23456L, "1.23456e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.234567L, "1.234567e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.2345678L, "1.2345678e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.23456789L, "1.23456789e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.234567890L, "1.23456789e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.2345678901L, "1.2345678901e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.23456789012L, "1.23456789012e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.234567890123L, "1.234567890123e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.2345678901234L, "1.2345678901234e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.23456789012345L, "1.23456789012345e+00", boost::charconv::chars_format::scientific);
    test_printf_fallback(1.234567890123456L, "1.234567890123456e+00", boost::charconv::chars_format::scientific);

    test_printf_fallback(1.0L, "1e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.2L, "1.2e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.23L, "1.23e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.234L, "1.234e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.2345L, "1.2345e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.23456L, "1.23456e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.234567L, "1.234567e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.2345678L, "1.2345678e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.23456789L, "1.23456789e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.234567890L, "1.23456789e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.2345678901L, "1.2345678901e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.23456789012L, "1.23456789012e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.234567890123L, "1.234567890123e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.2345678901234L, "1.2345678901234e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.23456789012345L, "1.23456789012345e+00", boost::charconv::chars_format::fixed);
    test_printf_fallback(1.234567890123456L, "1.234567890123456e+00", boost::charconv::chars_format::fixed);

    // Regressions or numbers that take >64 bits to represent correctly
    spot_check(9007199254740991.0, "9.007199254740991e+15", boost::charconv::chars_format::scientific);
    spot_check(9007199254740992.0, "9.007199254740992e+15", boost::charconv::chars_format::scientific);
    spot_check(123456789012345683968.0, "1.2345678901234568e+20", boost::charconv::chars_format::scientific);
    spot_check(1.9430376160308388E16, "1.9430376160308388e+16", boost::charconv::chars_format::scientific);
    spot_check(-6.9741824662760956E19, "-6.9741824662760956e+19", boost::charconv::chars_format::scientific);
    spot_check(4.3816050601147837E18, "4.3816050601147837e+18", boost::charconv::chars_format::scientific);
    spot_check(1.8531501765868567E21, "1.8531501765868567e+21", boost::charconv::chars_format::scientific);
    spot_check(-3.347727380279489E33, "-3.347727380279489e+33", boost::charconv::chars_format::scientific);
    spot_check(9.409340012568248E18, "9.409340012568248e+18", boost::charconv::chars_format::scientific);
    spot_check(4.708356024711512E18, "4.708356024711512e+18", boost::charconv::chars_format::scientific);
    spot_check(9.0608011534336E15, "9.0608011534336e+15", boost::charconv::chars_format::scientific);
    spot_check(2.989102097996E-312, "2.989102097996e-312", boost::charconv::chars_format::scientific);
    spot_check(1.18575755E-316, "1.18575755e-316", boost::charconv::chars_format::scientific);
    spot_check(4.940656E-318, "4.940656e-318", boost::charconv::chars_format::scientific);

    // Every power
    spot_check(1.7e+308, "1.7e+308", boost::charconv::chars_format::scientific);
    spot_check(1.7e+307, "1.7e+307", boost::charconv::chars_format::scientific);
    spot_check(1.7e+306, "1.7e+306", boost::charconv::chars_format::scientific);
    spot_check(1.7e+305, "1.7e+305", boost::charconv::chars_format::scientific);
    spot_check(1.7e+304, "1.7e+304", boost::charconv::chars_format::scientific);
    spot_check(1.7e+303, "1.7e+303", boost::charconv::chars_format::scientific);
    spot_check(1.7e+302, "1.7e+302", boost::charconv::chars_format::scientific);
    spot_check(1.7e+301, "1.7e+301", boost::charconv::chars_format::scientific);
    spot_check(1.7e+300, "1.7e+300", boost::charconv::chars_format::scientific);
    spot_check(1.7e+299, "1.7e+299", boost::charconv::chars_format::scientific);
    spot_check(1.7e+298, "1.7e+298", boost::charconv::chars_format::scientific);
    spot_check(1.7e+297, "1.7e+297", boost::charconv::chars_format::scientific);
    spot_check(1.7e+296, "1.7e+296", boost::charconv::chars_format::scientific);
    spot_check(1.7e+295, "1.7e+295", boost::charconv::chars_format::scientific);
    spot_check(1.7e+294, "1.7e+294", boost::charconv::chars_format::scientific);
    spot_check(1.7e+293, "1.7e+293", boost::charconv::chars_format::scientific);
    spot_check(1.7e+292, "1.7e+292", boost::charconv::chars_format::scientific);
    spot_check(1.7e+291, "1.7e+291", boost::charconv::chars_format::scientific);
    spot_check(1.7e+290, "1.7e+290", boost::charconv::chars_format::scientific);
    spot_check(1.7e+289, "1.7e+289", boost::charconv::chars_format::scientific);
    spot_check(1.7e+288, "1.7e+288", boost::charconv::chars_format::scientific);
    spot_check(1.7e+287, "1.7e+287", boost::charconv::chars_format::scientific);
    spot_check(1.7e+286, "1.7e+286", boost::charconv::chars_format::scientific);
    spot_check(1.7e+285, "1.7e+285", boost::charconv::chars_format::scientific);
    spot_check(1.7e+284, "1.7e+284", boost::charconv::chars_format::scientific);
    spot_check(1.7e+283, "1.7e+283", boost::charconv::chars_format::scientific);
    spot_check(1.7e+282, "1.7e+282", boost::charconv::chars_format::scientific);
    spot_check(1.7e+281, "1.7e+281", boost::charconv::chars_format::scientific);
    spot_check(1.7e+280, "1.7e+280", boost::charconv::chars_format::scientific);
    spot_check(1.7e+279, "1.7e+279", boost::charconv::chars_format::scientific);
    spot_check(1.7e+278, "1.7e+278", boost::charconv::chars_format::scientific);
    spot_check(1.7e+277, "1.7e+277", boost::charconv::chars_format::scientific);
    spot_check(1.7e+276, "1.7e+276", boost::charconv::chars_format::scientific);
    spot_check(1.7e+275, "1.7e+275", boost::charconv::chars_format::scientific);
    spot_check(1.7e+274, "1.7e+274", boost::charconv::chars_format::scientific);
    spot_check(1.7e+273, "1.7e+273", boost::charconv::chars_format::scientific);
    spot_check(1.7e+272, "1.7e+272", boost::charconv::chars_format::scientific);
    spot_check(1.7e+271, "1.7e+271", boost::charconv::chars_format::scientific);
    spot_check(1.7e+270, "1.7e+270", boost::charconv::chars_format::scientific);
    spot_check(1.7e+269, "1.7e+269", boost::charconv::chars_format::scientific);
    spot_check(1.7e+268, "1.7e+268", boost::charconv::chars_format::scientific);
    spot_check(1.7e+267, "1.7e+267", boost::charconv::chars_format::scientific);
    spot_check(1.7e+266, "1.7e+266", boost::charconv::chars_format::scientific);
    spot_check(1.7e+265, "1.7e+265", boost::charconv::chars_format::scientific);
    spot_check(1.7e+264, "1.7e+264", boost::charconv::chars_format::scientific);
    spot_check(1.7e+263, "1.7e+263", boost::charconv::chars_format::scientific);
    spot_check(1.7e+262, "1.7e+262", boost::charconv::chars_format::scientific);
    spot_check(1.7e+261, "1.7e+261", boost::charconv::chars_format::scientific);
    spot_check(1.7e+260, "1.7e+260", boost::charconv::chars_format::scientific);
    spot_check(1.7e+259, "1.7e+259", boost::charconv::chars_format::scientific);
    spot_check(1.7e+258, "1.7e+258", boost::charconv::chars_format::scientific);
    spot_check(1.7e+257, "1.7e+257", boost::charconv::chars_format::scientific);
    spot_check(1.7e+256, "1.7e+256", boost::charconv::chars_format::scientific);
    spot_check(1.7e+255, "1.7e+255", boost::charconv::chars_format::scientific);
    spot_check(1.7e+254, "1.7e+254", boost::charconv::chars_format::scientific);
    spot_check(1.7e+253, "1.7e+253", boost::charconv::chars_format::scientific);
    spot_check(1.7e+252, "1.7e+252", boost::charconv::chars_format::scientific);
    spot_check(1.7e+251, "1.7e+251", boost::charconv::chars_format::scientific);
    spot_check(1.7e+250, "1.7e+250", boost::charconv::chars_format::scientific);
    spot_check(1.7e+249, "1.7e+249", boost::charconv::chars_format::scientific);
    spot_check(1.7e+248, "1.7e+248", boost::charconv::chars_format::scientific);
    spot_check(1.7e+247, "1.7e+247", boost::charconv::chars_format::scientific);
    spot_check(1.7e+246, "1.7e+246", boost::charconv::chars_format::scientific);
    spot_check(1.7e+245, "1.7e+245", boost::charconv::chars_format::scientific);
    spot_check(1.7e+244, "1.7e+244", boost::charconv::chars_format::scientific);
    spot_check(1.7e+243, "1.7e+243", boost::charconv::chars_format::scientific);
    spot_check(1.7e+242, "1.7e+242", boost::charconv::chars_format::scientific);
    spot_check(1.7e+241, "1.7e+241", boost::charconv::chars_format::scientific);
    spot_check(1.7e+240, "1.7e+240", boost::charconv::chars_format::scientific);
    spot_check(1.7e+239, "1.7e+239", boost::charconv::chars_format::scientific);
    spot_check(1.7e+238, "1.7e+238", boost::charconv::chars_format::scientific);
    spot_check(1.7e+237, "1.7e+237", boost::charconv::chars_format::scientific);
    spot_check(1.7e+236, "1.7e+236", boost::charconv::chars_format::scientific);
    spot_check(1.7e+235, "1.7e+235", boost::charconv::chars_format::scientific);
    spot_check(1.7e+234, "1.7e+234", boost::charconv::chars_format::scientific);
    spot_check(1.7e+233, "1.7e+233", boost::charconv::chars_format::scientific);
    spot_check(1.7e+232, "1.7e+232", boost::charconv::chars_format::scientific);
    spot_check(1.7e+231, "1.7e+231", boost::charconv::chars_format::scientific);
    spot_check(1.7e+230, "1.7e+230", boost::charconv::chars_format::scientific);
    spot_check(1.7e+229, "1.7e+229", boost::charconv::chars_format::scientific);
    spot_check(1.7e+228, "1.7e+228", boost::charconv::chars_format::scientific);
    spot_check(1.7e+227, "1.7e+227", boost::charconv::chars_format::scientific);
    spot_check(1.7e+226, "1.7e+226", boost::charconv::chars_format::scientific);
    spot_check(1.7e+225, "1.7e+225", boost::charconv::chars_format::scientific);
    spot_check(1.7e+224, "1.7e+224", boost::charconv::chars_format::scientific);
    spot_check(1.7e+223, "1.7e+223", boost::charconv::chars_format::scientific);
    spot_check(1.7e+222, "1.7e+222", boost::charconv::chars_format::scientific);
    spot_check(1.7e+221, "1.7e+221", boost::charconv::chars_format::scientific);
    spot_check(1.7e+220, "1.7e+220", boost::charconv::chars_format::scientific);
    spot_check(1.7e+219, "1.7e+219", boost::charconv::chars_format::scientific);
    spot_check(1.7e+218, "1.7e+218", boost::charconv::chars_format::scientific);
    spot_check(1.7e+217, "1.7e+217", boost::charconv::chars_format::scientific);
    spot_check(1.7e+216, "1.7e+216", boost::charconv::chars_format::scientific);
    spot_check(1.7e+215, "1.7e+215", boost::charconv::chars_format::scientific);
    spot_check(1.7e+214, "1.7e+214", boost::charconv::chars_format::scientific);
    spot_check(1.7e+213, "1.7e+213", boost::charconv::chars_format::scientific);
    spot_check(1.7e+212, "1.7e+212", boost::charconv::chars_format::scientific);
    spot_check(1.7e+211, "1.7e+211", boost::charconv::chars_format::scientific);
    spot_check(1.7e+210, "1.7e+210", boost::charconv::chars_format::scientific);
    spot_check(1.7e+209, "1.7e+209", boost::charconv::chars_format::scientific);
    spot_check(1.7e+208, "1.7e+208", boost::charconv::chars_format::scientific);
    spot_check(1.7e+207, "1.7e+207", boost::charconv::chars_format::scientific);
    spot_check(1.7e+206, "1.7e+206", boost::charconv::chars_format::scientific);
    spot_check(1.7e+205, "1.7e+205", boost::charconv::chars_format::scientific);
    spot_check(1.7e+204, "1.7e+204", boost::charconv::chars_format::scientific);
    spot_check(1.7e+203, "1.7e+203", boost::charconv::chars_format::scientific);
    spot_check(1.7e+202, "1.7e+202", boost::charconv::chars_format::scientific);
    spot_check(1.7e+201, "1.7e+201", boost::charconv::chars_format::scientific);
    spot_check(1.7e+200, "1.7e+200", boost::charconv::chars_format::scientific);

    spot_check(1.7e+199, "1.7e+199", boost::charconv::chars_format::scientific);
    spot_check(1.7e+198, "1.7e+198", boost::charconv::chars_format::scientific);
    spot_check(1.7e+197, "1.7e+197", boost::charconv::chars_format::scientific);
    spot_check(1.7e+196, "1.7e+196", boost::charconv::chars_format::scientific);
    spot_check(1.7e+195, "1.7e+195", boost::charconv::chars_format::scientific);
    spot_check(1.7e+194, "1.7e+194", boost::charconv::chars_format::scientific);
    spot_check(1.7e+193, "1.7e+193", boost::charconv::chars_format::scientific);
    spot_check(1.7e+192, "1.7e+192", boost::charconv::chars_format::scientific);
    spot_check(1.7e+191, "1.7e+191", boost::charconv::chars_format::scientific);
    spot_check(1.7e+190, "1.7e+190", boost::charconv::chars_format::scientific);
    spot_check(1.7e+189, "1.7e+189", boost::charconv::chars_format::scientific);
    spot_check(1.7e+188, "1.7e+188", boost::charconv::chars_format::scientific);
    spot_check(1.7e+187, "1.7e+187", boost::charconv::chars_format::scientific);
    spot_check(1.7e+186, "1.7e+186", boost::charconv::chars_format::scientific);
    spot_check(1.7e+185, "1.7e+185", boost::charconv::chars_format::scientific);
    spot_check(1.7e+184, "1.7e+184", boost::charconv::chars_format::scientific);
    spot_check(1.7e+183, "1.7e+183", boost::charconv::chars_format::scientific);
    spot_check(1.7e+182, "1.7e+182", boost::charconv::chars_format::scientific);
    spot_check(1.7e+181, "1.7e+181", boost::charconv::chars_format::scientific);
    spot_check(1.7e+180, "1.7e+180", boost::charconv::chars_format::scientific);
    spot_check(1.7e+179, "1.7e+179", boost::charconv::chars_format::scientific);
    spot_check(1.7e+178, "1.7e+178", boost::charconv::chars_format::scientific);
    spot_check(1.7e+177, "1.7e+177", boost::charconv::chars_format::scientific);
    spot_check(1.7e+176, "1.7e+176", boost::charconv::chars_format::scientific);
    spot_check(1.7e+175, "1.7e+175", boost::charconv::chars_format::scientific);
    spot_check(1.7e+174, "1.7e+174", boost::charconv::chars_format::scientific);
    spot_check(1.7e+173, "1.7e+173", boost::charconv::chars_format::scientific);
    spot_check(1.7e+172, "1.7e+172", boost::charconv::chars_format::scientific);
    spot_check(1.7e+171, "1.7e+171", boost::charconv::chars_format::scientific);
    spot_check(1.7e+170, "1.7e+170", boost::charconv::chars_format::scientific);
    spot_check(1.7e+169, "1.7e+169", boost::charconv::chars_format::scientific);
    spot_check(1.7e+168, "1.7e+168", boost::charconv::chars_format::scientific);
    spot_check(1.7e+167, "1.7e+167", boost::charconv::chars_format::scientific);
    spot_check(1.7e+166, "1.7e+166", boost::charconv::chars_format::scientific);
    spot_check(1.7e+165, "1.7e+165", boost::charconv::chars_format::scientific);
    spot_check(1.7e+164, "1.7e+164", boost::charconv::chars_format::scientific);
    spot_check(1.7e+163, "1.7e+163", boost::charconv::chars_format::scientific);
    spot_check(1.7e+162, "1.7e+162", boost::charconv::chars_format::scientific);
    spot_check(1.7e+161, "1.7e+161", boost::charconv::chars_format::scientific);
    spot_check(1.7e+160, "1.7e+160", boost::charconv::chars_format::scientific);
    spot_check(1.7e+159, "1.7e+159", boost::charconv::chars_format::scientific);
    spot_check(1.7e+158, "1.7e+158", boost::charconv::chars_format::scientific);
    spot_check(1.7e+157, "1.7e+157", boost::charconv::chars_format::scientific);
    spot_check(1.7e+156, "1.7e+156", boost::charconv::chars_format::scientific);
    spot_check(1.7e+155, "1.7e+155", boost::charconv::chars_format::scientific);
    spot_check(1.7e+154, "1.7e+154", boost::charconv::chars_format::scientific);
    spot_check(1.7e+153, "1.7e+153", boost::charconv::chars_format::scientific);
    spot_check(1.7e+152, "1.7e+152", boost::charconv::chars_format::scientific);
    spot_check(1.7e+151, "1.7e+151", boost::charconv::chars_format::scientific);
    spot_check(1.7e+150, "1.7e+150", boost::charconv::chars_format::scientific);
    spot_check(1.7e+149, "1.7e+149", boost::charconv::chars_format::scientific);
    spot_check(1.7e+148, "1.7e+148", boost::charconv::chars_format::scientific);
    spot_check(1.7e+147, "1.7e+147", boost::charconv::chars_format::scientific);
    spot_check(1.7e+146, "1.7e+146", boost::charconv::chars_format::scientific);
    spot_check(1.7e+145, "1.7e+145", boost::charconv::chars_format::scientific);
    spot_check(1.7e+144, "1.7e+144", boost::charconv::chars_format::scientific);
    spot_check(1.7e+143, "1.7e+143", boost::charconv::chars_format::scientific);
    spot_check(1.7e+142, "1.7e+142", boost::charconv::chars_format::scientific);
    spot_check(1.7e+141, "1.7e+141", boost::charconv::chars_format::scientific);
    spot_check(1.7e+140, "1.7e+140", boost::charconv::chars_format::scientific);
    spot_check(1.7e+139, "1.7e+139", boost::charconv::chars_format::scientific);
    spot_check(1.7e+138, "1.7e+138", boost::charconv::chars_format::scientific);
    spot_check(1.7e+137, "1.7e+137", boost::charconv::chars_format::scientific);
    spot_check(1.7e+136, "1.7e+136", boost::charconv::chars_format::scientific);
    spot_check(1.7e+135, "1.7e+135", boost::charconv::chars_format::scientific);
    spot_check(1.7e+134, "1.7e+134", boost::charconv::chars_format::scientific);
    spot_check(1.7e+133, "1.7e+133", boost::charconv::chars_format::scientific);
    spot_check(1.7e+132, "1.7e+132", boost::charconv::chars_format::scientific);
    spot_check(1.7e+131, "1.7e+131", boost::charconv::chars_format::scientific);
    spot_check(1.7e+130, "1.7e+130", boost::charconv::chars_format::scientific);
    spot_check(1.7e+129, "1.7e+129", boost::charconv::chars_format::scientific);
    spot_check(1.7e+128, "1.7e+128", boost::charconv::chars_format::scientific);
    spot_check(1.7e+127, "1.7e+127", boost::charconv::chars_format::scientific);
    spot_check(1.7e+126, "1.7e+126", boost::charconv::chars_format::scientific);
    spot_check(1.7e+125, "1.7e+125", boost::charconv::chars_format::scientific);
    spot_check(1.7e+124, "1.7e+124", boost::charconv::chars_format::scientific);
    spot_check(1.7e+123, "1.7e+123", boost::charconv::chars_format::scientific);
    spot_check(1.7e+122, "1.7e+122", boost::charconv::chars_format::scientific);
    spot_check(1.7e+121, "1.7e+121", boost::charconv::chars_format::scientific);
    spot_check(1.7e+120, "1.7e+120", boost::charconv::chars_format::scientific);
    spot_check(1.7e+119, "1.7e+119", boost::charconv::chars_format::scientific);
    spot_check(1.7e+118, "1.7e+118", boost::charconv::chars_format::scientific);
    spot_check(1.7e+117, "1.7e+117", boost::charconv::chars_format::scientific);
    spot_check(1.7e+116, "1.7e+116", boost::charconv::chars_format::scientific);
    spot_check(1.7e+115, "1.7e+115", boost::charconv::chars_format::scientific);
    spot_check(1.7e+114, "1.7e+114", boost::charconv::chars_format::scientific);
    spot_check(1.7e+113, "1.7e+113", boost::charconv::chars_format::scientific);
    spot_check(1.7e+112, "1.7e+112", boost::charconv::chars_format::scientific);
    spot_check(1.7e+111, "1.7e+111", boost::charconv::chars_format::scientific);
    spot_check(1.7e+110, "1.7e+110", boost::charconv::chars_format::scientific);
    spot_check(1.7e+109, "1.7e+109", boost::charconv::chars_format::scientific);
    spot_check(1.7e+108, "1.7e+108", boost::charconv::chars_format::scientific);
    spot_check(1.7e+107, "1.7e+107", boost::charconv::chars_format::scientific);
    spot_check(1.7e+106, "1.7e+106", boost::charconv::chars_format::scientific);
    spot_check(1.7e+105, "1.7e+105", boost::charconv::chars_format::scientific);
    spot_check(1.7e+104, "1.7e+104", boost::charconv::chars_format::scientific);
    spot_check(1.7e+103, "1.7e+103", boost::charconv::chars_format::scientific);
    spot_check(1.7e+102, "1.7e+102", boost::charconv::chars_format::scientific);
    spot_check(1.7e+101, "1.7e+101", boost::charconv::chars_format::scientific);
    spot_check(1.7e+100, "1.7e+100", boost::charconv::chars_format::scientific);

    spot_check(1.7e+99, "1.7e+99", boost::charconv::chars_format::scientific);
    spot_check(1.7e+98, "1.7e+98", boost::charconv::chars_format::scientific);
    spot_check(1.7e+97, "1.7e+97", boost::charconv::chars_format::scientific);
    spot_check(1.7e+96, "1.7e+96", boost::charconv::chars_format::scientific);
    spot_check(1.7e+95, "1.7e+95", boost::charconv::chars_format::scientific);
    spot_check(1.7e+94, "1.7e+94", boost::charconv::chars_format::scientific);
    spot_check(1.7e+93, "1.7e+93", boost::charconv::chars_format::scientific);
    spot_check(1.7e+92, "1.7e+92", boost::charconv::chars_format::scientific);
    spot_check(1.7e+91, "1.7e+91", boost::charconv::chars_format::scientific);
    spot_check(1.7e+90, "1.7e+90", boost::charconv::chars_format::scientific);
    spot_check(1.7e+89, "1.7e+89", boost::charconv::chars_format::scientific);
    spot_check(1.7e+88, "1.7e+88", boost::charconv::chars_format::scientific);
    spot_check(1.7e+87, "1.7e+87", boost::charconv::chars_format::scientific);
    spot_check(1.7e+86, "1.7e+86", boost::charconv::chars_format::scientific);
    spot_check(1.7e+85, "1.7e+85", boost::charconv::chars_format::scientific);
    spot_check(1.7e+84, "1.7e+84", boost::charconv::chars_format::scientific);
    spot_check(1.7e+83, "1.7e+83", boost::charconv::chars_format::scientific);
    spot_check(1.7e+82, "1.7e+82", boost::charconv::chars_format::scientific);
    spot_check(1.7e+81, "1.7e+81", boost::charconv::chars_format::scientific);
    spot_check(1.7e+80, "1.7e+80", boost::charconv::chars_format::scientific);
    spot_check(1.7e+79, "1.7e+79", boost::charconv::chars_format::scientific);
    spot_check(1.7e+78, "1.7e+78", boost::charconv::chars_format::scientific);
    spot_check(1.7e+77, "1.7e+77", boost::charconv::chars_format::scientific);
    spot_check(1.7e+76, "1.7e+76", boost::charconv::chars_format::scientific);
    spot_check(1.7e+75, "1.7e+75", boost::charconv::chars_format::scientific);
    spot_check(1.7e+74, "1.7e+74", boost::charconv::chars_format::scientific);
    spot_check(1.7e+73, "1.7e+73", boost::charconv::chars_format::scientific);
    spot_check(1.7e+72, "1.7e+72", boost::charconv::chars_format::scientific);
    spot_check(1.7e+71, "1.7e+71", boost::charconv::chars_format::scientific);
    spot_check(1.7e+70, "1.7e+70", boost::charconv::chars_format::scientific);
    spot_check(1.7e+69, "1.7e+69", boost::charconv::chars_format::scientific);
    spot_check(1.7e+68, "1.7e+68", boost::charconv::chars_format::scientific);
    spot_check(1.7e+67, "1.7e+67", boost::charconv::chars_format::scientific);
    spot_check(1.7e+66, "1.7e+66", boost::charconv::chars_format::scientific);
    spot_check(1.7e+65, "1.7e+65", boost::charconv::chars_format::scientific);
    spot_check(1.7e+64, "1.7e+64", boost::charconv::chars_format::scientific);
    spot_check(1.7e+63, "1.7e+63", boost::charconv::chars_format::scientific);
    spot_check(1.7e+62, "1.7e+62", boost::charconv::chars_format::scientific);
    spot_check(1.7e+61, "1.7e+61", boost::charconv::chars_format::scientific);
    spot_check(1.7e+60, "1.7e+60", boost::charconv::chars_format::scientific);
    spot_check(1.7e+59, "1.7e+59", boost::charconv::chars_format::scientific);
    spot_check(1.7e+58, "1.7e+58", boost::charconv::chars_format::scientific);
    spot_check(1.7e+57, "1.7e+57", boost::charconv::chars_format::scientific);
    spot_check(1.7e+56, "1.7e+56", boost::charconv::chars_format::scientific);
    spot_check(1.7e+55, "1.7e+55", boost::charconv::chars_format::scientific);
    spot_check(1.7e+54, "1.7e+54", boost::charconv::chars_format::scientific);
    spot_check(1.7e+53, "1.7e+53", boost::charconv::chars_format::scientific);
    spot_check(1.7e+52, "1.7e+52", boost::charconv::chars_format::scientific);
    spot_check(1.7e+51, "1.7e+51", boost::charconv::chars_format::scientific);
    spot_check(1.7e+50, "1.7e+50", boost::charconv::chars_format::scientific);
    spot_check(1.7e+49, "1.7e+49", boost::charconv::chars_format::scientific);
    spot_check(1.7e+48, "1.7e+48", boost::charconv::chars_format::scientific);
    spot_check(1.7e+47, "1.7e+47", boost::charconv::chars_format::scientific);
    spot_check(1.7e+46, "1.7e+46", boost::charconv::chars_format::scientific);
    spot_check(1.7e+45, "1.7e+45", boost::charconv::chars_format::scientific);
    spot_check(1.7e+44, "1.7e+44", boost::charconv::chars_format::scientific);
    spot_check(1.7e+43, "1.7e+43", boost::charconv::chars_format::scientific);
    spot_check(1.7e+42, "1.7e+42", boost::charconv::chars_format::scientific);
    spot_check(1.7e+41, "1.7e+41", boost::charconv::chars_format::scientific);
    spot_check(1.7e+40, "1.7e+40", boost::charconv::chars_format::scientific);
    spot_check(1.7e+39, "1.7e+39", boost::charconv::chars_format::scientific);
    spot_check(1.7e+38, "1.7e+38", boost::charconv::chars_format::scientific);
    spot_check(1.7e+37, "1.7e+37", boost::charconv::chars_format::scientific);
    spot_check(1.7e+36, "1.7e+36", boost::charconv::chars_format::scientific);
    spot_check(1.7e+35, "1.7e+35", boost::charconv::chars_format::scientific);
    spot_check(1.7e+34, "1.7e+34", boost::charconv::chars_format::scientific);
    spot_check(1.7e+33, "1.7e+33", boost::charconv::chars_format::scientific);
    spot_check(1.7e+32, "1.7e+32", boost::charconv::chars_format::scientific);
    spot_check(1.7e+31, "1.7e+31", boost::charconv::chars_format::scientific);
    spot_check(1.7e+30, "1.7e+30", boost::charconv::chars_format::scientific);
    spot_check(1.7e+29, "1.7e+29", boost::charconv::chars_format::scientific);
    spot_check(1.7e+28, "1.7e+28", boost::charconv::chars_format::scientific);
    spot_check(1.7e+27, "1.7e+27", boost::charconv::chars_format::scientific);
    spot_check(1.7e+26, "1.7e+26", boost::charconv::chars_format::scientific);
    spot_check(1.7e+25, "1.7e+25", boost::charconv::chars_format::scientific);
    spot_check(1.7e+24, "1.7e+24", boost::charconv::chars_format::scientific);
    spot_check(1.7e+23, "1.7e+23", boost::charconv::chars_format::scientific);
    spot_check(1.7e+22, "1.7e+22", boost::charconv::chars_format::scientific);
    spot_check(1.7e+21, "1.7e+21", boost::charconv::chars_format::scientific);
    spot_check(1.7e+20, "1.7e+20", boost::charconv::chars_format::scientific);
    spot_check(1.7e+19, "1.7e+19", boost::charconv::chars_format::scientific);
    spot_check(1.7e+18, "1.7e+18", boost::charconv::chars_format::scientific);
    spot_check(1.7e+17, "1.7e+17", boost::charconv::chars_format::scientific);
    spot_check(1.7e+16, "1.7e+16", boost::charconv::chars_format::scientific);
    spot_check(1.7e+15, "1.7e+15", boost::charconv::chars_format::scientific);
    spot_check(1.7e+14, "1.7e+14", boost::charconv::chars_format::scientific);
    spot_check(1.7e+13, "1.7e+13", boost::charconv::chars_format::scientific);
    spot_check(1.7e+12, "1.7e+12", boost::charconv::chars_format::scientific);
    spot_check(1.7e+11, "1.7e+11", boost::charconv::chars_format::scientific);
    spot_check(1.7e+10, "1.7e+10", boost::charconv::chars_format::scientific);
    spot_check(1.7e+09, "1.7e+09", boost::charconv::chars_format::scientific);
    spot_check(1.7e+08, "1.7e+08", boost::charconv::chars_format::scientific);
    spot_check(1.7e+07, "1.7e+07", boost::charconv::chars_format::scientific);
    spot_check(1.7e+06, "1.7e+06", boost::charconv::chars_format::scientific);
    spot_check(1.7e+05, "1.7e+05", boost::charconv::chars_format::scientific);
    spot_check(1.7e+04, "1.7e+04", boost::charconv::chars_format::scientific);
    spot_check(1.7e+03, "1.7e+03", boost::charconv::chars_format::scientific);
    spot_check(1.7e+02, "1.7e+02", boost::charconv::chars_format::scientific);
    spot_check(1.7e+01, "1.7e+01", boost::charconv::chars_format::scientific);
    spot_check(1.7e+00, "1.7e+00", boost::charconv::chars_format::scientific);

    spot_check(1.7e-308, "1.7e-308", boost::charconv::chars_format::scientific);
    spot_check(1.7e-307, "1.7e-307", boost::charconv::chars_format::scientific);
    spot_check(1.7e-306, "1.7e-306", boost::charconv::chars_format::scientific);
    spot_check(1.7e-305, "1.7e-305", boost::charconv::chars_format::scientific);
    spot_check(1.7e-304, "1.7e-304", boost::charconv::chars_format::scientific);
    spot_check(1.7e-303, "1.7e-303", boost::charconv::chars_format::scientific);
    spot_check(1.7e-302, "1.7e-302", boost::charconv::chars_format::scientific);
    spot_check(1.7e-301, "1.7e-301", boost::charconv::chars_format::scientific);
    spot_check(1.7e-300, "1.7e-300", boost::charconv::chars_format::scientific);
    spot_check(1.7e-299, "1.7e-299", boost::charconv::chars_format::scientific);
    spot_check(1.7e-298, "1.7e-298", boost::charconv::chars_format::scientific);
    spot_check(1.7e-297, "1.7e-297", boost::charconv::chars_format::scientific);
    spot_check(1.7e-296, "1.7e-296", boost::charconv::chars_format::scientific);
    spot_check(1.7e-295, "1.7e-295", boost::charconv::chars_format::scientific);
    spot_check(1.7e-294, "1.7e-294", boost::charconv::chars_format::scientific);
    spot_check(1.7e-293, "1.7e-293", boost::charconv::chars_format::scientific);
    spot_check(1.7e-292, "1.7e-292", boost::charconv::chars_format::scientific);
    spot_check(1.7e-291, "1.7e-291", boost::charconv::chars_format::scientific);
    spot_check(1.7e-290, "1.7e-290", boost::charconv::chars_format::scientific);
    spot_check(1.7e-289, "1.7e-289", boost::charconv::chars_format::scientific);
    spot_check(1.7e-288, "1.7e-288", boost::charconv::chars_format::scientific);
    spot_check(1.7e-287, "1.7e-287", boost::charconv::chars_format::scientific);
    spot_check(1.7e-286, "1.7e-286", boost::charconv::chars_format::scientific);
    spot_check(1.7e-285, "1.7e-285", boost::charconv::chars_format::scientific);
    spot_check(1.7e-284, "1.7e-284", boost::charconv::chars_format::scientific);
    spot_check(1.7e-283, "1.7e-283", boost::charconv::chars_format::scientific);
    spot_check(1.7e-282, "1.7e-282", boost::charconv::chars_format::scientific);
    spot_check(1.7e-281, "1.7e-281", boost::charconv::chars_format::scientific);
    spot_check(1.7e-280, "1.7e-280", boost::charconv::chars_format::scientific);
    spot_check(1.7e-279, "1.7e-279", boost::charconv::chars_format::scientific);
    spot_check(1.7e-278, "1.7e-278", boost::charconv::chars_format::scientific);
    spot_check(1.7e-277, "1.7e-277", boost::charconv::chars_format::scientific);
    spot_check(1.7e-276, "1.7e-276", boost::charconv::chars_format::scientific);
    spot_check(1.7e-275, "1.7e-275", boost::charconv::chars_format::scientific);
    spot_check(1.7e-274, "1.7e-274", boost::charconv::chars_format::scientific);
    spot_check(1.7e-273, "1.7e-273", boost::charconv::chars_format::scientific);
    spot_check(1.7e-272, "1.7e-272", boost::charconv::chars_format::scientific);
    spot_check(1.7e-271, "1.7e-271", boost::charconv::chars_format::scientific);
    spot_check(1.7e-270, "1.7e-270", boost::charconv::chars_format::scientific);
    spot_check(1.7e-269, "1.7e-269", boost::charconv::chars_format::scientific);
    spot_check(1.7e-268, "1.7e-268", boost::charconv::chars_format::scientific);
    spot_check(1.7e-267, "1.7e-267", boost::charconv::chars_format::scientific);
    spot_check(1.7e-266, "1.7e-266", boost::charconv::chars_format::scientific);
    spot_check(1.7e-265, "1.7e-265", boost::charconv::chars_format::scientific);
    spot_check(1.7e-264, "1.7e-264", boost::charconv::chars_format::scientific);
    spot_check(1.7e-263, "1.7e-263", boost::charconv::chars_format::scientific);
    spot_check(1.7e-262, "1.7e-262", boost::charconv::chars_format::scientific);
    spot_check(1.7e-261, "1.7e-261", boost::charconv::chars_format::scientific);
    spot_check(1.7e-260, "1.7e-260", boost::charconv::chars_format::scientific);
    spot_check(1.7e-259, "1.7e-259", boost::charconv::chars_format::scientific);
    spot_check(1.7e-258, "1.7e-258", boost::charconv::chars_format::scientific);
    spot_check(1.7e-257, "1.7e-257", boost::charconv::chars_format::scientific);
    spot_check(1.7e-256, "1.7e-256", boost::charconv::chars_format::scientific);
    spot_check(1.7e-255, "1.7e-255", boost::charconv::chars_format::scientific);
    spot_check(1.7e-254, "1.7e-254", boost::charconv::chars_format::scientific);
    spot_check(1.7e-253, "1.7e-253", boost::charconv::chars_format::scientific);
    spot_check(1.7e-252, "1.7e-252", boost::charconv::chars_format::scientific);
    spot_check(1.7e-251, "1.7e-251", boost::charconv::chars_format::scientific);
    spot_check(1.7e-250, "1.7e-250", boost::charconv::chars_format::scientific);
    spot_check(1.7e-249, "1.7e-249", boost::charconv::chars_format::scientific);
    spot_check(1.7e-248, "1.7e-248", boost::charconv::chars_format::scientific);
    spot_check(1.7e-247, "1.7e-247", boost::charconv::chars_format::scientific);
    spot_check(1.7e-246, "1.7e-246", boost::charconv::chars_format::scientific);
    spot_check(1.7e-245, "1.7e-245", boost::charconv::chars_format::scientific);
    spot_check(1.7e-244, "1.7e-244", boost::charconv::chars_format::scientific);
    spot_check(1.7e-243, "1.7e-243", boost::charconv::chars_format::scientific);
    spot_check(1.7e-242, "1.7e-242", boost::charconv::chars_format::scientific);
    spot_check(1.7e-241, "1.7e-241", boost::charconv::chars_format::scientific);
    spot_check(1.7e-240, "1.7e-240", boost::charconv::chars_format::scientific);
    spot_check(1.7e-239, "1.7e-239", boost::charconv::chars_format::scientific);
    spot_check(1.7e-238, "1.7e-238", boost::charconv::chars_format::scientific);
    spot_check(1.7e-237, "1.7e-237", boost::charconv::chars_format::scientific);
    spot_check(1.7e-236, "1.7e-236", boost::charconv::chars_format::scientific);
    spot_check(1.7e-235, "1.7e-235", boost::charconv::chars_format::scientific);
    spot_check(1.7e-234, "1.7e-234", boost::charconv::chars_format::scientific);
    spot_check(1.7e-233, "1.7e-233", boost::charconv::chars_format::scientific);
    spot_check(1.7e-232, "1.7e-232", boost::charconv::chars_format::scientific);
    spot_check(1.7e-231, "1.7e-231", boost::charconv::chars_format::scientific);
    spot_check(1.7e-230, "1.7e-230", boost::charconv::chars_format::scientific);
    spot_check(1.7e-229, "1.7e-229", boost::charconv::chars_format::scientific);
    spot_check(1.7e-228, "1.7e-228", boost::charconv::chars_format::scientific);
    spot_check(1.7e-227, "1.7e-227", boost::charconv::chars_format::scientific);
    spot_check(1.7e-226, "1.7e-226", boost::charconv::chars_format::scientific);
    spot_check(1.7e-225, "1.7e-225", boost::charconv::chars_format::scientific);
    spot_check(1.7e-224, "1.7e-224", boost::charconv::chars_format::scientific);
    spot_check(1.7e-223, "1.7e-223", boost::charconv::chars_format::scientific);
    spot_check(1.7e-222, "1.7e-222", boost::charconv::chars_format::scientific);
    spot_check(1.7e-221, "1.7e-221", boost::charconv::chars_format::scientific);
    spot_check(1.7e-220, "1.7e-220", boost::charconv::chars_format::scientific);
    spot_check(1.7e-219, "1.7e-219", boost::charconv::chars_format::scientific);
    spot_check(1.7e-218, "1.7e-218", boost::charconv::chars_format::scientific);
    spot_check(1.7e-217, "1.7e-217", boost::charconv::chars_format::scientific);
    spot_check(1.7e-216, "1.7e-216", boost::charconv::chars_format::scientific);
    spot_check(1.7e-215, "1.7e-215", boost::charconv::chars_format::scientific);
    spot_check(1.7e-214, "1.7e-214", boost::charconv::chars_format::scientific);
    spot_check(1.7e-213, "1.7e-213", boost::charconv::chars_format::scientific);
    spot_check(1.7e-212, "1.7e-212", boost::charconv::chars_format::scientific);
    spot_check(1.7e-211, "1.7e-211", boost::charconv::chars_format::scientific);
    spot_check(1.7e-210, "1.7e-210", boost::charconv::chars_format::scientific);
    spot_check(1.7e-209, "1.7e-209", boost::charconv::chars_format::scientific);
    spot_check(1.7e-208, "1.7e-208", boost::charconv::chars_format::scientific);
    spot_check(1.7e-207, "1.7e-207", boost::charconv::chars_format::scientific);
    spot_check(1.7e-206, "1.7e-206", boost::charconv::chars_format::scientific);
    spot_check(1.7e-205, "1.7e-205", boost::charconv::chars_format::scientific);
    spot_check(1.7e-204, "1.7e-204", boost::charconv::chars_format::scientific);
    spot_check(1.7e-203, "1.7e-203", boost::charconv::chars_format::scientific);
    spot_check(1.7e-202, "1.7e-202", boost::charconv::chars_format::scientific);
    spot_check(1.7e-201, "1.7e-201", boost::charconv::chars_format::scientific);
    spot_check(1.7e-200, "1.7e-200", boost::charconv::chars_format::scientific);

    spot_check(1.7e-199, "1.7e-199", boost::charconv::chars_format::scientific);
    spot_check(1.7e-198, "1.7e-198", boost::charconv::chars_format::scientific);
    spot_check(1.7e-197, "1.7e-197", boost::charconv::chars_format::scientific);
    spot_check(1.7e-196, "1.7e-196", boost::charconv::chars_format::scientific);
    spot_check(1.7e-195, "1.7e-195", boost::charconv::chars_format::scientific);
    spot_check(1.7e-194, "1.7e-194", boost::charconv::chars_format::scientific);
    spot_check(1.7e-193, "1.7e-193", boost::charconv::chars_format::scientific);
    spot_check(1.7e-192, "1.7e-192", boost::charconv::chars_format::scientific);
    spot_check(1.7e-191, "1.7e-191", boost::charconv::chars_format::scientific);
    spot_check(1.7e-190, "1.7e-190", boost::charconv::chars_format::scientific);
    spot_check(1.7e-189, "1.7e-189", boost::charconv::chars_format::scientific);
    spot_check(1.7e-188, "1.7e-188", boost::charconv::chars_format::scientific);
    spot_check(1.7e-187, "1.7e-187", boost::charconv::chars_format::scientific);
    spot_check(1.7e-186, "1.7e-186", boost::charconv::chars_format::scientific);
    spot_check(1.7e-185, "1.7e-185", boost::charconv::chars_format::scientific);
    spot_check(1.7e-184, "1.7e-184", boost::charconv::chars_format::scientific);
    spot_check(1.7e-183, "1.7e-183", boost::charconv::chars_format::scientific);
    spot_check(1.7e-182, "1.7e-182", boost::charconv::chars_format::scientific);
    spot_check(1.7e-181, "1.7e-181", boost::charconv::chars_format::scientific);
    spot_check(1.7e-180, "1.7e-180", boost::charconv::chars_format::scientific);
    spot_check(1.7e-179, "1.7e-179", boost::charconv::chars_format::scientific);
    spot_check(1.7e-178, "1.7e-178", boost::charconv::chars_format::scientific);
    spot_check(1.7e-177, "1.7e-177", boost::charconv::chars_format::scientific);
    spot_check(1.7e-176, "1.7e-176", boost::charconv::chars_format::scientific);
    spot_check(1.7e-175, "1.7e-175", boost::charconv::chars_format::scientific);
    spot_check(1.7e-174, "1.7e-174", boost::charconv::chars_format::scientific);
    spot_check(1.7e-173, "1.7e-173", boost::charconv::chars_format::scientific);
    spot_check(1.7e-172, "1.7e-172", boost::charconv::chars_format::scientific);
    spot_check(1.7e-171, "1.7e-171", boost::charconv::chars_format::scientific);
    spot_check(1.7e-170, "1.7e-170", boost::charconv::chars_format::scientific);
    spot_check(1.7e-169, "1.7e-169", boost::charconv::chars_format::scientific);
    spot_check(1.7e-168, "1.7e-168", boost::charconv::chars_format::scientific);
    spot_check(1.7e-167, "1.7e-167", boost::charconv::chars_format::scientific);
    spot_check(1.7e-166, "1.7e-166", boost::charconv::chars_format::scientific);
    spot_check(1.7e-165, "1.7e-165", boost::charconv::chars_format::scientific);
    spot_check(1.7e-164, "1.7e-164", boost::charconv::chars_format::scientific);
    spot_check(1.7e-163, "1.7e-163", boost::charconv::chars_format::scientific);
    spot_check(1.7e-162, "1.7e-162", boost::charconv::chars_format::scientific);
    spot_check(1.7e-161, "1.7e-161", boost::charconv::chars_format::scientific);
    spot_check(1.7e-160, "1.7e-160", boost::charconv::chars_format::scientific);
    spot_check(1.7e-159, "1.7e-159", boost::charconv::chars_format::scientific);
    spot_check(1.7e-158, "1.7e-158", boost::charconv::chars_format::scientific);
    spot_check(1.7e-157, "1.7e-157", boost::charconv::chars_format::scientific);
    spot_check(1.7e-156, "1.7e-156", boost::charconv::chars_format::scientific);
    spot_check(1.7e-155, "1.7e-155", boost::charconv::chars_format::scientific);
    spot_check(1.7e-154, "1.7e-154", boost::charconv::chars_format::scientific);
    spot_check(1.7e-153, "1.7e-153", boost::charconv::chars_format::scientific);
    spot_check(1.7e-152, "1.7e-152", boost::charconv::chars_format::scientific);
    spot_check(1.7e-151, "1.7e-151", boost::charconv::chars_format::scientific);
    spot_check(1.7e-150, "1.7e-150", boost::charconv::chars_format::scientific);
    spot_check(1.7e-149, "1.7e-149", boost::charconv::chars_format::scientific);
    spot_check(1.7e-148, "1.7e-148", boost::charconv::chars_format::scientific);
    spot_check(1.7e-147, "1.7e-147", boost::charconv::chars_format::scientific);
    spot_check(1.7e-146, "1.7e-146", boost::charconv::chars_format::scientific);
    spot_check(1.7e-145, "1.7e-145", boost::charconv::chars_format::scientific);
    spot_check(1.7e-144, "1.7e-144", boost::charconv::chars_format::scientific);
    spot_check(1.7e-143, "1.7e-143", boost::charconv::chars_format::scientific);
    spot_check(1.7e-142, "1.7e-142", boost::charconv::chars_format::scientific);
    spot_check(1.7e-141, "1.7e-141", boost::charconv::chars_format::scientific);
    spot_check(1.7e-140, "1.7e-140", boost::charconv::chars_format::scientific);
    spot_check(1.7e-139, "1.7e-139", boost::charconv::chars_format::scientific);
    spot_check(1.7e-138, "1.7e-138", boost::charconv::chars_format::scientific);
    spot_check(1.7e-137, "1.7e-137", boost::charconv::chars_format::scientific);
    spot_check(1.7e-136, "1.7e-136", boost::charconv::chars_format::scientific);
    spot_check(1.7e-135, "1.7e-135", boost::charconv::chars_format::scientific);
    spot_check(1.7e-134, "1.7e-134", boost::charconv::chars_format::scientific);
    spot_check(1.7e-133, "1.7e-133", boost::charconv::chars_format::scientific);
    spot_check(1.7e-132, "1.7e-132", boost::charconv::chars_format::scientific);
    spot_check(1.7e-131, "1.7e-131", boost::charconv::chars_format::scientific);
    spot_check(1.7e-130, "1.7e-130", boost::charconv::chars_format::scientific);
    spot_check(1.7e-129, "1.7e-129", boost::charconv::chars_format::scientific);
    spot_check(1.7e-128, "1.7e-128", boost::charconv::chars_format::scientific);
    spot_check(1.7e-127, "1.7e-127", boost::charconv::chars_format::scientific);
    spot_check(1.7e-126, "1.7e-126", boost::charconv::chars_format::scientific);
    spot_check(1.7e-125, "1.7e-125", boost::charconv::chars_format::scientific);
    spot_check(1.7e-124, "1.7e-124", boost::charconv::chars_format::scientific);
    spot_check(1.7e-123, "1.7e-123", boost::charconv::chars_format::scientific);
    spot_check(1.7e-122, "1.7e-122", boost::charconv::chars_format::scientific);
    spot_check(1.7e-121, "1.7e-121", boost::charconv::chars_format::scientific);
    spot_check(1.7e-120, "1.7e-120", boost::charconv::chars_format::scientific);
    spot_check(1.7e-119, "1.7e-119", boost::charconv::chars_format::scientific);
    spot_check(1.7e-118, "1.7e-118", boost::charconv::chars_format::scientific);
    spot_check(1.7e-117, "1.7e-117", boost::charconv::chars_format::scientific);
    spot_check(1.7e-116, "1.7e-116", boost::charconv::chars_format::scientific);
    spot_check(1.7e-115, "1.7e-115", boost::charconv::chars_format::scientific);
    spot_check(1.7e-114, "1.7e-114", boost::charconv::chars_format::scientific);
    spot_check(1.7e-113, "1.7e-113", boost::charconv::chars_format::scientific);
    spot_check(1.7e-112, "1.7e-112", boost::charconv::chars_format::scientific);
    spot_check(1.7e-111, "1.7e-111", boost::charconv::chars_format::scientific);
    spot_check(1.7e-110, "1.7e-110", boost::charconv::chars_format::scientific);
    spot_check(1.7e-109, "1.7e-109", boost::charconv::chars_format::scientific);
    spot_check(1.7e-108, "1.7e-108", boost::charconv::chars_format::scientific);
    spot_check(1.7e-107, "1.7e-107", boost::charconv::chars_format::scientific);
    spot_check(1.7e-106, "1.7e-106", boost::charconv::chars_format::scientific);
    spot_check(1.7e-105, "1.7e-105", boost::charconv::chars_format::scientific);
    spot_check(1.7e-104, "1.7e-104", boost::charconv::chars_format::scientific);
    spot_check(1.7e-103, "1.7e-103", boost::charconv::chars_format::scientific);
    spot_check(1.7e-102, "1.7e-102", boost::charconv::chars_format::scientific);
    spot_check(1.7e-101, "1.7e-101", boost::charconv::chars_format::scientific);
    spot_check(1.7e-100, "1.7e-100", boost::charconv::chars_format::scientific);

    spot_check(1.7e-99, "1.7e-99", boost::charconv::chars_format::scientific);
    spot_check(1.7e-98, "1.7e-98", boost::charconv::chars_format::scientific);
    spot_check(1.7e-97, "1.7e-97", boost::charconv::chars_format::scientific);
    spot_check(1.7e-96, "1.7e-96", boost::charconv::chars_format::scientific);
    spot_check(1.7e-95, "1.7e-95", boost::charconv::chars_format::scientific);
    spot_check(1.7e-94, "1.7e-94", boost::charconv::chars_format::scientific);
    spot_check(1.7e-93, "1.7e-93", boost::charconv::chars_format::scientific);
    spot_check(1.7e-92, "1.7e-92", boost::charconv::chars_format::scientific);
    spot_check(1.7e-91, "1.7e-91", boost::charconv::chars_format::scientific);
    spot_check(1.7e-90, "1.7e-90", boost::charconv::chars_format::scientific);
    spot_check(1.7e-89, "1.7e-89", boost::charconv::chars_format::scientific);
    spot_check(1.7e-88, "1.7e-88", boost::charconv::chars_format::scientific);
    spot_check(1.7e-87, "1.7e-87", boost::charconv::chars_format::scientific);
    spot_check(1.7e-86, "1.7e-86", boost::charconv::chars_format::scientific);
    spot_check(1.7e-85, "1.7e-85", boost::charconv::chars_format::scientific);
    spot_check(1.7e-84, "1.7e-84", boost::charconv::chars_format::scientific);
    spot_check(1.7e-83, "1.7e-83", boost::charconv::chars_format::scientific);
    spot_check(1.7e-82, "1.7e-82", boost::charconv::chars_format::scientific);
    spot_check(1.7e-81, "1.7e-81", boost::charconv::chars_format::scientific);
    spot_check(1.7e-80, "1.7e-80", boost::charconv::chars_format::scientific);
    spot_check(1.7e-79, "1.7e-79", boost::charconv::chars_format::scientific);
    spot_check(1.7e-78, "1.7e-78", boost::charconv::chars_format::scientific);
    spot_check(1.7e-77, "1.7e-77", boost::charconv::chars_format::scientific);
    spot_check(1.7e-76, "1.7e-76", boost::charconv::chars_format::scientific);
    spot_check(1.7e-75, "1.7e-75", boost::charconv::chars_format::scientific);
    spot_check(1.7e-74, "1.7e-74", boost::charconv::chars_format::scientific);
    spot_check(1.7e-73, "1.7e-73", boost::charconv::chars_format::scientific);
    spot_check(1.7e-72, "1.7e-72", boost::charconv::chars_format::scientific);
    spot_check(1.7e-71, "1.7e-71", boost::charconv::chars_format::scientific);
    spot_check(1.7e-70, "1.7e-70", boost::charconv::chars_format::scientific);
    spot_check(1.7e-69, "1.7e-69", boost::charconv::chars_format::scientific);
    spot_check(1.7e-68, "1.7e-68", boost::charconv::chars_format::scientific);
    spot_check(1.7e-67, "1.7e-67", boost::charconv::chars_format::scientific);
    spot_check(1.7e-66, "1.7e-66", boost::charconv::chars_format::scientific);
    spot_check(1.7e-65, "1.7e-65", boost::charconv::chars_format::scientific);
    spot_check(1.7e-64, "1.7e-64", boost::charconv::chars_format::scientific);
    spot_check(1.7e-63, "1.7e-63", boost::charconv::chars_format::scientific);
    spot_check(1.7e-62, "1.7e-62", boost::charconv::chars_format::scientific);
    spot_check(1.7e-61, "1.7e-61", boost::charconv::chars_format::scientific);
    spot_check(1.7e-60, "1.7e-60", boost::charconv::chars_format::scientific);
    spot_check(1.7e-59, "1.7e-59", boost::charconv::chars_format::scientific);
    spot_check(1.7e-58, "1.7e-58", boost::charconv::chars_format::scientific);
    spot_check(1.7e-57, "1.7e-57", boost::charconv::chars_format::scientific);
    spot_check(1.7e-56, "1.7e-56", boost::charconv::chars_format::scientific);
    spot_check(1.7e-55, "1.7e-55", boost::charconv::chars_format::scientific);
    spot_check(1.7e-54, "1.7e-54", boost::charconv::chars_format::scientific);
    spot_check(1.7e-53, "1.7e-53", boost::charconv::chars_format::scientific);
    spot_check(1.7e-52, "1.7e-52", boost::charconv::chars_format::scientific);
    spot_check(1.7e-51, "1.7e-51", boost::charconv::chars_format::scientific);
    spot_check(1.7e-50, "1.7e-50", boost::charconv::chars_format::scientific);
    spot_check(1.7e-49, "1.7e-49", boost::charconv::chars_format::scientific);
    spot_check(1.7e-48, "1.7e-48", boost::charconv::chars_format::scientific);
    spot_check(1.7e-47, "1.7e-47", boost::charconv::chars_format::scientific);
    spot_check(1.7e-46, "1.7e-46", boost::charconv::chars_format::scientific);
    spot_check(1.7e-45, "1.7e-45", boost::charconv::chars_format::scientific);
    spot_check(1.7e-44, "1.7e-44", boost::charconv::chars_format::scientific);
    spot_check(1.7e-43, "1.7e-43", boost::charconv::chars_format::scientific);
    spot_check(1.7e-42, "1.7e-42", boost::charconv::chars_format::scientific);
    spot_check(1.7e-41, "1.7e-41", boost::charconv::chars_format::scientific);
    spot_check(1.7e-40, "1.7e-40", boost::charconv::chars_format::scientific);
    spot_check(1.7e-39, "1.7e-39", boost::charconv::chars_format::scientific);
    spot_check(1.7e-38, "1.7e-38", boost::charconv::chars_format::scientific);
    spot_check(1.7e-37, "1.7e-37", boost::charconv::chars_format::scientific);
    spot_check(1.7e-36, "1.7e-36", boost::charconv::chars_format::scientific);
    spot_check(1.7e-35, "1.7e-35", boost::charconv::chars_format::scientific);
    spot_check(1.7e-34, "1.7e-34", boost::charconv::chars_format::scientific);
    spot_check(1.7e-33, "1.7e-33", boost::charconv::chars_format::scientific);
    spot_check(1.7e-32, "1.7e-32", boost::charconv::chars_format::scientific);
    spot_check(1.7e-31, "1.7e-31", boost::charconv::chars_format::scientific);
    spot_check(1.7e-30, "1.7e-30", boost::charconv::chars_format::scientific);
    spot_check(1.7e-29, "1.7e-29", boost::charconv::chars_format::scientific);
    spot_check(1.7e-28, "1.7e-28", boost::charconv::chars_format::scientific);
    spot_check(1.7e-27, "1.7e-27", boost::charconv::chars_format::scientific);
    spot_check(1.7e-26, "1.7e-26", boost::charconv::chars_format::scientific);
    spot_check(1.7e-25, "1.7e-25", boost::charconv::chars_format::scientific);
    spot_check(1.7e-24, "1.7e-24", boost::charconv::chars_format::scientific);
    spot_check(1.7e-23, "1.7e-23", boost::charconv::chars_format::scientific);
    spot_check(1.7e-22, "1.7e-22", boost::charconv::chars_format::scientific);
    spot_check(1.7e-21, "1.7e-21", boost::charconv::chars_format::scientific);
    spot_check(1.7e-20, "1.7e-20", boost::charconv::chars_format::scientific);
    spot_check(1.7e-19, "1.7e-19", boost::charconv::chars_format::scientific);
    spot_check(1.7e-18, "1.7e-18", boost::charconv::chars_format::scientific);
    spot_check(1.7e-17, "1.7e-17", boost::charconv::chars_format::scientific);
    spot_check(1.7e-16, "1.7e-16", boost::charconv::chars_format::scientific);
    spot_check(1.7e-15, "1.7e-15", boost::charconv::chars_format::scientific);
    spot_check(1.7e-14, "1.7e-14", boost::charconv::chars_format::scientific);
    spot_check(1.7e-13, "1.7e-13", boost::charconv::chars_format::scientific);
    spot_check(1.7e-12, "1.7e-12", boost::charconv::chars_format::scientific);
    spot_check(1.7e-11, "1.7e-11", boost::charconv::chars_format::scientific);
    spot_check(1.7e-10, "1.7e-10", boost::charconv::chars_format::scientific);
    spot_check(1.7e-09, "1.7e-09", boost::charconv::chars_format::scientific);
    spot_check(1.7e-08, "1.7e-08", boost::charconv::chars_format::scientific);
    spot_check(1.7e-07, "1.7e-07", boost::charconv::chars_format::scientific);
    spot_check(1.7e-06, "1.7e-06", boost::charconv::chars_format::scientific);
    spot_check(1.7e-05, "1.7e-05", boost::charconv::chars_format::scientific);
    spot_check(1.7e-04, "1.7e-04", boost::charconv::chars_format::scientific);
    spot_check(1.7e-03, "1.7e-03", boost::charconv::chars_format::scientific);
    spot_check(1.7e-02, "1.7e-02", boost::charconv::chars_format::scientific);
    spot_check(1.7e-01, "1.7e-01", boost::charconv::chars_format::scientific);
    spot_check(1.7e-00, "1.7e+00", boost::charconv::chars_format::scientific);

    spot_check(0.0, "0", boost::charconv::chars_format::fixed);
    spot_check(-0.0, "-0", boost::charconv::chars_format::fixed);
    spot_check(0.0, "0.0000000000", boost::charconv::chars_format::fixed, 10);
    spot_check(-0.0, "-0.0000000000", boost::charconv::chars_format::fixed, 10);

    #ifdef BOOST_CHARCONV_HAS_FLOAT16
    {
        constexpr int N = 1024;
        boost::detail::splitmix64 rng(42);

        std::float16_t const small_q = std::pow(1.0F16, -16.0F16);
        constexpr std::uint64_t divisor = UINT64_MAX / 16;

        for( int i = 0; i < N; ++i )
        {
            std::float16_t w0 = static_cast<std::float16_t>( rng() / divisor );
            test_to_chars_hex_16(w0);

            std::float16_t w1 = static_cast<std::float16_t>( rng() / divisor ) * small_q;
            test_to_chars_hex_16(w1);

            std::float16_t w2 = (std::numeric_limits<std::float16_t>::max)() / static_cast<std::float16_t>( rng() / divisor );
            test_to_chars_hex_16(w2);

            std::float16_t w3 = (std::numeric_limits<std::float16_t>::min)() * static_cast<std::float16_t>( rng() / divisor );
            test_to_chars_hex_16(w3);
        }
    }
    #endif

    #ifdef BOOST_CHARCONV_HAS_BRAINFLOAT16
    {
        constexpr int N = 1024;
        boost::detail::splitmix64 rng(42);

        std::bfloat16_t const small_q = std::pow(1.0BF16, -16.0BF16);

        for( int i = 0; i < N; ++i )
        {
            std::bfloat16_t w0 = static_cast<std::bfloat16_t>( rng() );
            test_to_chars_hex_16(w0);

            std::bfloat16_t w1 = static_cast<std::bfloat16_t>( rng() ) * small_q;
            test_to_chars_hex_16(w1);

            std::bfloat16_t w2 = (std::numeric_limits<std::bfloat16_t>::max)() / static_cast<std::bfloat16_t>( rng() );
            test_to_chars_hex_16(w2);

            std::bfloat16_t w3 = (std::numeric_limits<std::bfloat16_t>::min)() * static_cast<std::bfloat16_t>( rng() );
            test_to_chars_hex_16(w3);
        }
    }
    #endif

    return boost::report_errors();
}
