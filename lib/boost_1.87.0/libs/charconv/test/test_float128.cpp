// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/charconv/detail/config.hpp>

#if defined(BOOST_CHARCONV_HAS_QUADMATH) && defined(BOOST_HAS_INT128)

#include <ostream>

#ifdef BOOST_CHARCONV_HAS_STDFLOAT128
#include <charconv>

std::ostream& operator<<( std::ostream& os, __float128 v )
{
    char buffer[ 256 ] {};
    std::to_chars(buffer, buffer + sizeof(buffer), static_cast<std::float128_t>(v));
    os << buffer;
    return os;
}

std::ostream& operator<<( std::ostream& os, std::float128_t v)
{
    char buffer [ 256 ] {};
    std::to_chars(buffer, buffer + sizeof(buffer), v);
    os << buffer;
    return os;
}

#else

std::ostream& operator<<( std::ostream& os, __float128 v )
{
    char buffer[ 256 ] {};
    quadmath_snprintf(buffer, sizeof(buffer), "%Qg", v);
    os << buffer;
    return os;
}

#endif // BOOST_CHARCONV_HAS_STDFLOAT128

static char* mini_to_chars( char (&buffer)[ 64 ], boost::uint128_type v )
{
    char* p = buffer + 64;
    *--p = '\0';

    do
    {
        *--p = "0123456789"[ v % 10 ];
        v /= 10;
    }
    while ( v != 0 );

    return p;
}

std::ostream& operator<<( std::ostream& os, boost::uint128_type v )
{
    char buffer[ 64 ];

    os << mini_to_chars( buffer, v );
    return os;
}

std::ostream& operator<<( std::ostream& os, boost::int128_type v )
{
    char buffer[ 64 ];
    char* p;

    if( v >= 0 )
    {
        p = mini_to_chars( buffer, static_cast<boost::uint128_type>(v) );
    }
    else
    {
        p = mini_to_chars( buffer, -static_cast<boost::uint128_type>(v) );
        *--p = '-';
    }

    os << p;
    return os;
}


#include <boost/charconv.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/core/detail/splitmix64.hpp>
#include <boost/charconv/detail/generate_nan.hpp>
#include <boost/charconv/detail/issignaling.hpp>
#include <limits>
#include <iostream>
#include <iomanip>
#include <string>
#include <random>

constexpr int N = 1024;
static boost::detail::splitmix64 rng;

template <typename T>
void test_signaling_nan()
{
    BOOST_IF_CONSTEXPR (std::is_same<T, __float128>::value)
    {
        #if BOOST_CHARCONV_HAS_BUILTIN(__builtin_nansq)
        BOOST_TEST(boost::charconv::detail::issignaling(__builtin_nansq("")));
        BOOST_TEST(boost::charconv::detail::issignaling(-__builtin_nansq("")));
        #endif
    }
    else
    {
        #ifdef BOOST_CHARCONV_HAS_STDFLOAT128
        BOOST_TEST(boost::charconv::detail::issignaling(std::numeric_limits<T>::signaling_NaN()));
        BOOST_TEST(boost::charconv::detail::issignaling(-std::numeric_limits<T>::signaling_NaN()));
        #endif
    }

    BOOST_TEST(!(boost::charconv::detail::issignaling)(std::numeric_limits<T>::quiet_NaN()));
    BOOST_TEST(!(boost::charconv::detail::issignaling)(std::numeric_limits<T>::infinity()));
    BOOST_TEST(!(boost::charconv::detail::issignaling)(-std::numeric_limits<T>::quiet_NaN()));
    BOOST_TEST(!(boost::charconv::detail::issignaling)(-std::numeric_limits<T>::infinity()));
}

template <typename T>
boost::int128_type float_distance(T a, T b)
{
    boost::int128_type ai;
    boost::int128_type bi;

    std::memcpy(&ai, &a, sizeof(__float128));
    std::memcpy(&bi, &b, sizeof(__float128));

    boost::int128_type result = bi - ai;

    if (ai < 0 || bi < 0)
    {
        result = -result;
    }

    return result;
}

template <typename T>
void overflow_spot_value(const std::string& buffer, T expected_value, boost::charconv::chars_format fmt = boost::charconv::chars_format::general)
{
    auto v = static_cast<T>(42.Q);
    auto r = boost::charconv::from_chars_erange(buffer.c_str(), buffer.c_str() + std::strlen(buffer.c_str()), v, fmt);

    if (!(BOOST_TEST(v == expected_value) && BOOST_TEST(r.ec == std::errc::result_out_of_range)))
    {
        std::cerr << "Test failure for: " << buffer << " got: " << v << std::endl;
    }
}

template <typename T>
void test_issue_37()
{
    overflow_spot_value("1e99999", HUGE_VALQ);
    overflow_spot_value("-1e99999",-HUGE_VALQ);
    overflow_spot_value("1.0e+99999", HUGE_VALQ);
    overflow_spot_value("-1.0e+99999", -HUGE_VALQ);

    overflow_spot_value("1e-99999", static_cast<T>(0.0Q));
    overflow_spot_value("-1.0e-99999", static_cast<T>(-0.0Q));
}

template <typename T>
void test_roundtrip( T value )
{
    char buffer[ 256 ];

    auto r = boost::charconv::to_chars( buffer, buffer + sizeof( buffer ), value );

    BOOST_TEST( r.ec == std::errc() );

    T v2 = 0;
    auto r2 = boost::charconv::from_chars( buffer, r.ptr, v2 );

    if( BOOST_TEST( r2.ec == std::errc() ) && BOOST_TEST( std::abs(float_distance(v2, value)) <= 1 ) && BOOST_TEST( r2.ptr == r.ptr) )
    {
    }
    else
    {
        std::cerr << std::setprecision(35)
                  << "     Value: " << value
                  << "\n  To chars: " << std::string( buffer, r.ptr )
                  << "\nFrom chars: " << v2
                  << "\nULP distance: " << float_distance(v2, value) << std::endl;
    }
}

template <typename T>
const char* fmt_from_type (T)
{
    return "%Qg";
}

template <typename T>
const char* fmt_from_type_fixed (T)
{
    return "%.0f";
}

template <typename T>
const char* fmt_from_type_scientific (T)
{
    return "%.35Qe";
}

template <typename T>
const char* fmt_from_type_hex (T)
{
    return "%Qa";
}

inline int print(__float128 value, char* buffer, size_t buffer_size, const char* fmt)
{
    return quadmath_snprintf(buffer, buffer_size, fmt, value);
}

#ifdef BOOST_CHARCONV_HAS_STDFLOAT128
// Has no overload of strtod/sprintf etc so cast to __float128
// See: https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2022/p1467r9.html#library
inline int print(std::float128_t value, char* buffer, size_t buffer_size, const char* fmt)
{
    return print(static_cast<__float128>(value), buffer, buffer_size, fmt);
}

#endif

template <typename T>
void test_sprintf_float( T value, boost::charconv::chars_format fmt = boost::charconv::chars_format::scientific )
{
    char buffer [ 256 ] {};

    if (fmt == boost::charconv::chars_format::fixed && (value > 1e100L || value < 1e-100L))
    {
        // Avoid failures from overflow
        return;
    }

    const auto r = boost::charconv::to_chars( buffer, buffer + sizeof(buffer), value, fmt );
    BOOST_TEST( r.ec == std::errc() );

    char buffer2 [ 256 ] {};

    const char* sprintf_fmt;
    const char* error_format;
    switch (fmt)
    {
        case boost::charconv::chars_format::general:
            sprintf_fmt = fmt_from_type(value);
            error_format = "General";
            break;
        case boost::charconv::chars_format::scientific:
            sprintf_fmt = fmt_from_type_scientific(value);
            error_format = "Scientific";
            break;
        case boost::charconv::chars_format::fixed:
            sprintf_fmt = fmt_from_type_fixed(value);
            error_format = "Fixed";
            break;
        case boost::charconv::chars_format::hex:
            sprintf_fmt = fmt_from_type_hex(value);
            error_format = "Hex";
            break;
    }

    print( value, buffer2, sizeof(buffer2), sprintf_fmt );

    // Remove trailing zeros from printf (if applicable)
    std::string printf_string {buffer2};
    if (fmt == boost::charconv::chars_format::scientific)
    {
        std::size_t found_trailing_0 = printf_string.find_first_of('e');
        if (found_trailing_0 != std::string::npos)
        {
            --found_trailing_0;
            while (printf_string[found_trailing_0] == '0')
            {
                printf_string.erase(found_trailing_0, 1);
                --found_trailing_0;
            }
        }
    }
    else if (fmt == boost::charconv::chars_format::hex)
    {
        printf_string.erase(0, 2); // Remove 0x that printf appends
    }

    // Same issues that arise in to_chars_snprintf.cpp so abort if in range
    //
    //     Value: 3.350549627872214798203501062446534e-4913
    //  To chars: 3.350549627872214798203501062446534e-4913
    //  Snprintf: 3.35055e-4913
    //
    //     Value: 6.8220421318020332664117517756596e+4913
    //  To chars: 6.8220421318020332664117517756596e+4913
    //  Snprintf: 6.82204e+4913
    //
    //     Value: 1.0600979293241972185e-109
    //  To chars: 1.0600979293241972185e-109
    //  Snprintf: 1.0601e-109
    //
    if ((value > static_cast<T>(1e15Q) && value < static_cast<T>(1e20Q)) ||
        (value > static_cast<T>(1e4912Q) || value < static_cast<T>(1e-4912Q)) ||
        (value > static_cast<T>(1e-115Q) && value < static_cast<T>(2e-109Q)))
    {
        return;
    }

    if( BOOST_TEST_EQ( std::string( buffer, r.ptr ), printf_string ) )
    {
    }
    else
    {
        std::cerr << std::setprecision(35)
                  << "     Value: " << value
                  << "\n  To chars: " << std::string( buffer, r.ptr )
                  << "\n  Snprintf: " << printf_string
                  << "\n    Format: " << error_format << std::endl;
    }
}

template <typename T>
void test_roundtrip_bv()
{
    test_roundtrip( static_cast<T>(FLT128_MIN) );
    test_roundtrip( static_cast<T>(-FLT128_MIN) );
    test_roundtrip( static_cast<T>(FLT128_MAX) );
    test_roundtrip( static_cast<T>(-FLT128_MAX) );
}

#ifdef BOOST_CHARCONV_HAS_STDFLOAT128
template <typename T>
void test_spot(T val, boost::charconv::chars_format fmt = boost::charconv::chars_format::general, int precision = -1)
{
    if (fmt == boost::charconv::chars_format::fixed && (val > 1e100 || val < 1e-100))
    {
        // Avoid failres from overflow
        return;
    }

    std::chars_format stl_fmt;
    switch (fmt)
    {
        case boost::charconv::chars_format::general:
            stl_fmt = std::chars_format::general;
            break;
        case boost::charconv::chars_format::fixed:
            stl_fmt = std::chars_format::fixed;
            break;
        case boost::charconv::chars_format::scientific:
            stl_fmt = std::chars_format::scientific;
            break;
        case boost::charconv::chars_format::hex:
            stl_fmt = std::chars_format::hex;
            break;
        default:
            BOOST_UNREACHABLE_RETURN(fmt);
            break;
    }

    char buffer_boost[256];
    char buffer_stl[256];

    boost::charconv::to_chars_result r_boost;
    std::to_chars_result r_stl;

    if (precision == -1)
    {
        r_boost = boost::charconv::to_chars(buffer_boost, buffer_boost + sizeof(buffer_boost), val, fmt);
        r_stl = std::to_chars(buffer_stl, buffer_stl + sizeof(buffer_stl), val, stl_fmt);
    }
    else
    {
        r_boost = boost::charconv::to_chars(buffer_boost, buffer_boost + sizeof(buffer_boost), val, fmt, precision);
        r_stl = std::to_chars(buffer_stl, buffer_stl + sizeof(buffer_stl), val, stl_fmt, precision);
    }

    BOOST_TEST(r_boost.ec == std::errc());
    if (r_stl.ec != std::errc())
    {
        // STL failed
        return;
    }

    const std::ptrdiff_t diff_boost = r_boost.ptr - buffer_boost;
    const std::ptrdiff_t diff_stl = r_stl.ptr - buffer_stl;
    const auto boost_str = std::string(buffer_boost, r_boost.ptr);
    const auto stl_str = std::string(buffer_stl, r_stl.ptr);

    // Region of divergence between our results and that of the STL
    // Value: 13501897678889699601
    // Boost: 13501897678889699601
    //   STL: 1.3501897678889699601e+19
    if (val > static_cast<T>(1e15) && val < static_cast<T>(1e20))
    {
        return;
    }

    if (!(BOOST_TEST_CSTR_EQ(boost_str.c_str(), stl_str.c_str()) && BOOST_TEST_EQ(diff_boost, diff_stl)))
    {
        std::cerr << std::setprecision(35)
                  << "Value: " << val
                  << "\nBoost: " << boost_str.c_str()
                  << "\n  STL: " << stl_str.c_str() << std::endl;
    }
}

template <typename T>
void random_test(boost::charconv::chars_format fmt = boost::charconv::chars_format::general)
{
    std::mt19937_64 gen(42);
    std::uniform_real_distribution<T> dist(0, FLT128_MAX);

    for (int i = 0; i < N/2; ++i)
    {
        test_spot<T>(dist(gen), fmt);
    }

    // Test small values
    std::uniform_real_distribution<T> small_dist(0, 1);
    for (int i = 0; i < N/2; ++i)
    {
        test_spot<T>(small_dist(gen), fmt);
    }
}

boost::int128_type float_total = 0;
boost::uint128_type abs_float_total = 0;

template <typename T>
void charconv_roundtrip(T val, boost::charconv::chars_format fmt = boost::charconv::chars_format::general, int precision = -1)
{
    if (fmt == boost::charconv::chars_format::fixed && (val > 1e100 || val < 1e-100))
    {
        // Avoid failres from overflow
        return;
    }

    std::chars_format stl_fmt;
    switch (fmt)
    {
        case boost::charconv::chars_format::general:
            stl_fmt = std::chars_format::general;
            break;
        case boost::charconv::chars_format::fixed:
            stl_fmt = std::chars_format::fixed;
            break;
        case boost::charconv::chars_format::scientific:
            stl_fmt = std::chars_format::scientific;
            break;
        case boost::charconv::chars_format::hex:
            stl_fmt = std::chars_format::hex;
            break;
        default:
            BOOST_UNREACHABLE_RETURN(fmt);
            break;
    }

    char buffer_boost[256];
    char buffer_stl[256];

    boost::charconv::to_chars_result r_boost;
    std::to_chars_result r_stl;

    if (precision == -1)
    {
        r_boost = boost::charconv::to_chars(buffer_boost, buffer_boost + sizeof(buffer_boost), val, fmt);
        r_stl = std::to_chars(buffer_stl, buffer_stl + sizeof(buffer_stl), val, stl_fmt);
    }
    else
    {
        r_boost = boost::charconv::to_chars(buffer_boost, buffer_boost + sizeof(buffer_boost), val, fmt, precision);
        r_stl = std::to_chars(buffer_stl, buffer_stl + sizeof(buffer_stl), val, stl_fmt, precision);
    }

    BOOST_TEST(r_boost.ec == std::errc());
    if (r_stl.ec != std::errc())
    {
        // STL failed
        return;
    }

    const std::ptrdiff_t diff_boost = r_boost.ptr - buffer_boost;
    const std::ptrdiff_t diff_stl = r_stl.ptr - buffer_stl;
    const auto boost_str = std::string(buffer_boost, r_boost.ptr);
    const auto stl_str = std::string(buffer_stl, r_stl.ptr);

    // Region of divergence between our results and that of the STL
    // Value: 13501897678889699601
    // Boost: 13501897678889699601
    //   STL: 1.3501897678889699601e+19
    if (val > static_cast<T>(1e15) && val < static_cast<T>(1e20))
    {
        return;
    }

    if (!(BOOST_TEST_CSTR_EQ(boost_str.c_str(), stl_str.c_str()) && BOOST_TEST_EQ(diff_boost, diff_stl)))
    {
        std::cerr << std::setprecision(35)
                  << "Value: " << val
                  << "\nBoost: " << boost_str.c_str()
                  << "\n  STL: " << stl_str.c_str() << std::endl;
    }

    T val_boost;
    T val_stl;
    const auto r_boost2 = boost::charconv::from_chars(buffer_boost, r_boost.ptr, val_boost, fmt);
    const auto r_stl2 = std::from_chars(buffer_stl, r_stl.ptr, val_stl, stl_fmt);

    BOOST_TEST(r_boost2.ec == std::errc());
    if (r_stl2.ec != std::errc())
    {
        // STL failed
        return;
    }

    if (BOOST_TEST( val_boost == val_stl ))
    {
    }
    else
    {
        const boost::int128_type dist = float_distance(val, val_boost);

        float_total += dist;
        abs_float_total += static_cast<boost::uint128_type>(dist < 0 ? -dist : dist);
        std::cerr << std::setprecision(35)
                  << "Value: " << val
                  << "\nBoost: " << val_boost
                  << "\n ULPs: " << dist
                  << "\n  STL: " << val_stl << std::endl;
    }
}

template <typename T>
void random_roundtrip(boost::charconv::chars_format fmt = boost::charconv::chars_format::general)
{
    std::mt19937_64 gen(42);
    std::uniform_real_distribution<T> dist(0, FLT128_MAX);

    for (int i = 0; i < N; ++i)
    {
        charconv_roundtrip<T>(dist(gen), fmt);
    }
}

#endif // BOOST_CHARCONV_HAS_STDFLOAT128

void spot_check_nan(const std::string& buffer, boost::charconv::chars_format fmt)
{
    __float128 v {};
    auto r = boost::charconv::from_chars(buffer.c_str(), buffer.c_str() + buffer.size(), v, fmt);
    if (!(BOOST_TEST(isnanq(v)) && BOOST_TEST(r)))
    {
        std::cerr << "Test failure for: " << buffer << " got: " << v << std::endl; // LCOV_EXCL_LINE
    }
}

void spot_check_inf(const std::string& buffer, boost::charconv::chars_format fmt)
{
    __float128 v {};
    auto r = boost::charconv::from_chars(buffer.c_str(), buffer.c_str() + buffer.size(), v, fmt);
    if (!(BOOST_TEST(isinfq(v)) && BOOST_TEST(r)))
    {
        std::cerr << "Test failure for: " << buffer << " got: " << v << std::endl; // LCOV_EXCL_LINE
    }
}

#if defined(__GNUC__) && __GNUC__ < 9 && __GNUC__ >= 5
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#endif

void test_nanq()
{
    BOOST_TEST(isnanq(boost::charconv::detail::nanq())) && BOOST_TEST(!issignaling(boost::charconv::detail::nanq()));
}

void test_nans()
{
    BOOST_TEST(isnanq(boost::charconv::detail::nans())) && BOOST_TEST(issignaling(boost::charconv::detail::nans()));
}

#if defined(__GNUC__) && __GNUC__ < 9 && __GNUC__ >= 5
#pragma GCC diagnostic pop
#endif

int main()
{
    #if BOOST_CHARCONV_LDBL_BITS == 128
    test_signaling_nan<long double>();

    // 128-bit long double
    {
        const long double q = powq( 1.0L, -128.0L );

        for( int i = 0; i < N; ++i )
        {
            long double w0 = static_cast<long double>( rng() ); // 0 .. 2^128
            test_roundtrip( w0 );
            test_sprintf_float( w0, boost::charconv::chars_format::general );
            test_sprintf_float( w0, boost::charconv::chars_format::scientific );
            test_sprintf_float( w0, boost::charconv::chars_format::fixed );
            test_sprintf_float( w0, boost::charconv::chars_format::hex );

            long double w1 = static_cast<long double>( rng() * q ); // 0.0 .. 1.0
            test_roundtrip( w1 );
            test_sprintf_float( w1, boost::charconv::chars_format::general );
            test_sprintf_float( w1, boost::charconv::chars_format::scientific );
            test_sprintf_float( w1, boost::charconv::chars_format::fixed );
            test_sprintf_float( w1, boost::charconv::chars_format::hex );

            long double w2 = FLT128_MAX / static_cast<long double>( rng() ); // large values
            test_roundtrip( w2 );
            test_sprintf_float( w2, boost::charconv::chars_format::general );
            test_sprintf_float( w2, boost::charconv::chars_format::scientific );
            test_sprintf_float( w2, boost::charconv::chars_format::fixed );
            test_sprintf_float( w2, boost::charconv::chars_format::hex );

            long double w3 = FLT128_MIN * static_cast<long double>( rng() ); // small values
            test_roundtrip( w3 );
            test_sprintf_float( w3, boost::charconv::chars_format::general );
            test_sprintf_float( w3, boost::charconv::chars_format::scientific );
            test_sprintf_float( w3, boost::charconv::chars_format::fixed );
            test_sprintf_float( w3, boost::charconv::chars_format::hex );
        }

        test_roundtrip_bv<__float128>();
    }
    #endif

    test_signaling_nan<__float128>();

    // __float128
    {
        const __float128 q = powq( 1.0Q, -128.0Q );

        for( int i = 0; i < N; ++i )
        {
            __float128 w0 = static_cast<__float128>( rng() ); // 0 .. 2^128
            test_roundtrip( w0 );
            test_sprintf_float( w0, boost::charconv::chars_format::general );
            test_sprintf_float( w0, boost::charconv::chars_format::scientific );
            test_sprintf_float( w0, boost::charconv::chars_format::fixed );
            test_sprintf_float( w0, boost::charconv::chars_format::hex );

            __float128 w1 = static_cast<__float128>( rng() * q ); // 0.0 .. 1.0
            test_roundtrip( w1 );
            test_sprintf_float( w1, boost::charconv::chars_format::general );
            test_sprintf_float( w1, boost::charconv::chars_format::scientific );
            test_sprintf_float( w1, boost::charconv::chars_format::fixed );
            test_sprintf_float( w1, boost::charconv::chars_format::hex );

            __float128 w2 = FLT128_MAX / static_cast<__float128>( rng() ); // large values
            test_roundtrip( w2 );
            test_sprintf_float( w2, boost::charconv::chars_format::general );
            test_sprintf_float( w2, boost::charconv::chars_format::scientific );
            test_sprintf_float( w2, boost::charconv::chars_format::fixed );
            test_sprintf_float( w2, boost::charconv::chars_format::hex );

            __float128 w3 = FLT128_MIN * static_cast<__float128>( rng() ); // small values
            test_roundtrip( w3 );
            test_sprintf_float( w3, boost::charconv::chars_format::general );
            test_sprintf_float( w3, boost::charconv::chars_format::scientific );
            test_sprintf_float( w3, boost::charconv::chars_format::fixed );
            test_sprintf_float( w3, boost::charconv::chars_format::hex );

            __float128 w5 = -static_cast<__float128>( rng() ); // -0 .. 2^128
            test_roundtrip( w5 );
            test_sprintf_float( w5, boost::charconv::chars_format::general );
            test_sprintf_float( w5, boost::charconv::chars_format::scientific );
            test_sprintf_float( w5, boost::charconv::chars_format::fixed );
            test_sprintf_float( w5, boost::charconv::chars_format::hex );
        }

        test_roundtrip_bv<__float128>();
    }
    
    #ifdef BOOST_CHARCONV_HAS_STDFLOAT128
    test_signaling_nan<std::float128_t>();

    // std::float128_t
    {
        const std::float128_t q = 1.0e-128F128;

        for( int i = 0; i < N; ++i )
        {
            std::float128_t w0 = static_cast<std::float128_t>( rng() ); // 0 .. 2^128
            test_roundtrip( w0 );
            charconv_roundtrip( w0 );
            test_sprintf_float( w0, boost::charconv::chars_format::general );
            test_sprintf_float( w0, boost::charconv::chars_format::scientific );
            test_sprintf_float( w0, boost::charconv::chars_format::fixed );
            test_sprintf_float( w0, boost::charconv::chars_format::hex );

            test_spot( w0, boost::charconv::chars_format::general );
            test_spot( w0, boost::charconv::chars_format::scientific );
            test_spot( w0, boost::charconv::chars_format::fixed );
            test_spot( w0, boost::charconv::chars_format::hex );

            test_spot( w0, boost::charconv::chars_format::general, 6 );
            test_spot( w0, boost::charconv::chars_format::scientific, 8 );

            std::float128_t w1 = static_cast<std::float128_t>( rng() * q ); // 0.0 .. 1.0
            test_roundtrip( w1 );
            charconv_roundtrip( w1 );
            test_sprintf_float( w1, boost::charconv::chars_format::general );
            test_sprintf_float( w1, boost::charconv::chars_format::scientific );
            test_sprintf_float( w1, boost::charconv::chars_format::fixed );
            test_sprintf_float( w1, boost::charconv::chars_format::hex );

            test_spot( w1, boost::charconv::chars_format::general );
            test_spot( w1, boost::charconv::chars_format::scientific );
            test_spot( w1, boost::charconv::chars_format::fixed );
            test_spot( w1, boost::charconv::chars_format::hex );

            test_spot( w1, boost::charconv::chars_format::general, 6 );
            test_spot( w1, boost::charconv::chars_format::scientific, 8 );

            // std::numeric_limits<std::float128_t> was not specialized until GCC-14
            // same with __float128

            std::float128_t w2 = static_cast<std::float128_t>(FLT128_MAX) / static_cast<std::float128_t>( rng() ); // large values
            test_roundtrip( w2 );
            charconv_roundtrip( w2 );
            test_sprintf_float( w2, boost::charconv::chars_format::general );
            test_sprintf_float( w2, boost::charconv::chars_format::scientific );
            test_sprintf_float( w2, boost::charconv::chars_format::fixed );
            test_sprintf_float( w2, boost::charconv::chars_format::hex );

            test_spot( w2, boost::charconv::chars_format::general );
            test_spot( w2, boost::charconv::chars_format::scientific );
            test_spot( w2, boost::charconv::chars_format::fixed );
            test_spot( w2, boost::charconv::chars_format::hex );

            test_spot( w2, boost::charconv::chars_format::general, 6 );
            test_spot( w2, boost::charconv::chars_format::scientific, 8 );

            std::float128_t w3 = static_cast<std::float128_t>(FLT128_MIN) * static_cast<std::float128_t>( rng() ); // small values
            test_roundtrip( w3 );
            charconv_roundtrip( w3 );
            test_sprintf_float( w3, boost::charconv::chars_format::general );
            test_sprintf_float( w3, boost::charconv::chars_format::scientific );
            test_sprintf_float( w3, boost::charconv::chars_format::fixed );
            test_sprintf_float( w3, boost::charconv::chars_format::hex );

            test_spot( w3, boost::charconv::chars_format::general );
            test_spot( w3, boost::charconv::chars_format::scientific );
            test_spot( w3, boost::charconv::chars_format::fixed );
            test_spot( w3, boost::charconv::chars_format::hex );

            test_spot( w3, boost::charconv::chars_format::general, 6 );
            test_spot( w3, boost::charconv::chars_format::scientific, 8 );
        }

        test_roundtrip_bv<std::float128_t>();
    }

    random_test<__float128>();
    random_test<std::float128_t>();
    random_roundtrip<std::float128_t>();

    test_issue_37<__float128>();
    test_issue_37<std::float128_t>();

    // Issue 64
    // Some of these are commented out because our answers differ from the STL
    //
    // Value: 0.001
    // Boost: 1e-03
    //   STL: 0.001
    //
    // Value: 1e-04
    // Boost: 1e-04
    //   STL: 0.0001

    test_spot<std::float128_t>(1e-01F128);
    test_spot<std::float128_t>(1e-02F128);
    // test_spot<std::float128_t>(1e-03F128);
    // test_spot<std::float128_t>(1e-04F128);
    test_spot<std::float128_t>(1.01e-01F128);
    test_spot<std::float128_t>(1.001e-01F128);
    test_spot<std::float128_t>(1.0001e-01F128);
    test_spot<std::float128_t>(1.00001e-01F128);
    test_spot<std::float128_t>(1.000001e-01F128);
    test_spot<std::float128_t>(1.0000001e-01F128);
    test_spot<std::float128_t>(1.01e-02F128);
    test_spot<std::float128_t>(1.001e-02F128);
    test_spot<std::float128_t>(1.0001e-02F128);
    test_spot<std::float128_t>(1.00001e-02F128);
    test_spot<std::float128_t>(1.000001e-02F128);
    test_spot<std::float128_t>(1.0000001e-02F128);
    test_spot<std::float128_t>(1.01e-03F128);
    test_spot<std::float128_t>(1.001e-03F128);
    test_spot<std::float128_t>(1.0001e-03F128);
    test_spot<std::float128_t>(1.00001e-03F128);
    test_spot<std::float128_t>(1.000001e-03F128);
    test_spot<std::float128_t>(1.0000001e-03F128);
    test_spot<std::float128_t>(1.01e-04F128);
    test_spot<std::float128_t>(1.001e-04F128);
    test_spot<std::float128_t>(1.0001e-04F128);
    test_spot<std::float128_t>(1.00001e-04F128);
    test_spot<std::float128_t>(1.000001e-04F128);
    test_spot<std::float128_t>(1.0000001e-04F128);
    test_spot<std::float128_t>(-3.589653987658756543653653365436e+04F128, boost::charconv::chars_format::hex);

    if (abs_float_total != 0)
    {
        std::cerr << std::setprecision(5)
                  << "\nAverage ULP distance: " << static_cast<double>(float_total) / static_cast<double>(N*4)
                  << "\nAbsolute ULP avererage: " << static_cast<double>(abs_float_total) / static_cast<double>(N*4)
                  << "\nTotal ULP distance: " << abs_float_total << std::endl;
    }

    #endif

    spot_check_nan("nan", boost::charconv::chars_format::general);
    spot_check_nan("-nan", boost::charconv::chars_format::general);
    spot_check_inf("inf", boost::charconv::chars_format::general);
    spot_check_inf("-inf", boost::charconv::chars_format::general);
    spot_check_nan("NAN", boost::charconv::chars_format::general);
    spot_check_nan("-NAN", boost::charconv::chars_format::general);
    spot_check_inf("INF", boost::charconv::chars_format::general);
    spot_check_inf("-INF", boost::charconv::chars_format::general);
    spot_check_nan("nan(snan)", boost::charconv::chars_format::general);
    spot_check_nan("-nan(snan)", boost::charconv::chars_format::general);

    test_nanq();
    #if defined(__GNUC__) && __GNUC__ >= 6
    test_nans();
    #endif

    return boost::report_errors();
}

#else

int main()
{
    return 0;
}

#endif
