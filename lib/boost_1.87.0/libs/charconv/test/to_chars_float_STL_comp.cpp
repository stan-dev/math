// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/config.hpp>

// https://en.cppreference.com/w/cpp/compiler_support/17
#if (defined(__GNUC__) && __GNUC__ >= 11) || \
    ((defined(__clang__) && __clang_major__ >= 14 && !defined(__APPLE__)) || (defined(__clang__) && defined(__APPLE__) && __clang_major__ >= 16)) || \
    (defined(_MSC_VER) && _MSC_VER >= 1924)

#include <boost/charconv.hpp>
#include <boost/core/lightweight_test.hpp>
#include <system_error>
#include <charconv>
#include <type_traits>
#include <limits>
#include <random>
#include <string>
#include <iomanip>
#include <iostream>
#include <cstring>
#include <cstdint>
#include <cerrno>

template <typename T>
void test_spot(T val, boost::charconv::chars_format fmt = boost::charconv::chars_format::general, int precision = -1)
{
    std::chars_format stl_fmt;
    std::string fmt_str;
    switch (fmt)
    {
    case boost::charconv::chars_format::general:
        stl_fmt = std::chars_format::general;
        fmt_str = "general";
        break;
    case boost::charconv::chars_format::fixed:
        stl_fmt = std::chars_format::fixed;
        fmt_str = "fixed";
        break;
    case boost::charconv::chars_format::scientific:
        stl_fmt = std::chars_format::scientific;
        fmt_str = "scientific";
        break;
    case boost::charconv::chars_format::hex:
        stl_fmt = std::chars_format::hex;
        fmt_str = "hex";
        break;
    // LCOV_EXCL_START
    default:
        BOOST_UNREACHABLE_RETURN(fmt);
        break;
    // LCOV_EXCL_STOP
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
        return; // LCOV_EXCL_LINE
    }

    const std::ptrdiff_t diff_boost = r_boost.ptr - buffer_boost;
    const std::ptrdiff_t diff_stl = r_stl.ptr - buffer_stl;
    *r_stl.ptr = '\0';
    *r_boost.ptr = '\0';
    constexpr T max_value = std::is_same<T, float>::value ? static_cast<T>(1e33F) : static_cast<T>(1e302);

    if (val > max_value)
    {
        return;
    }
    else if (!(BOOST_TEST_CSTR_EQ(buffer_boost, buffer_stl) && BOOST_TEST_EQ(diff_boost, diff_stl)))
    {
        // LCOV_EXCL_START
        std::cerr << std::setprecision(std::numeric_limits<T>::max_digits10 + 1)
                    << "Value: " << val
                    << "\nPrecision: " << precision
                    << "\nFormat: " << fmt_str
                    << "\nBoost: " << buffer_boost
                    << "\n  STL: " << buffer_stl << std::endl;
        // LCOV_EXCL_STOP
    }
}

template <typename T>
void random_test(boost::charconv::chars_format fmt = boost::charconv::chars_format::general, T min_val = 0, T max_val = std::numeric_limits<T>::max())
{   
    std::mt19937_64 gen(42);

    std::uniform_real_distribution<T> dist(min_val, max_val);

    for (int i = -1; i < std::numeric_limits<T>::digits10; ++i)
    {
        // TODO(mborland): fix integer part rounding
        if (fmt == boost::charconv::chars_format::hex && (i == 0 || i == 1 || i == 2))
        {
            continue;
        }

        for (std::size_t j = 0; j < 1000; ++j)
        {
            test_spot(dist(gen), fmt, i);
        }
    }

    // Test some extended precision arguments
    if (fmt != boost::charconv::chars_format::fixed || max_val < T(1e10))
    {
        for (int i = 40; i < 50; ++i)
        {
            for (std::size_t j = 0; j < 10; ++j)
            {
                test_spot(dist(gen), fmt, i);
            }
        }
    }
}

template <typename T>
void non_finite_test(boost::charconv::chars_format fmt = boost::charconv::chars_format::general)
{
    for (int i = -1; i <= 0; ++i)
    {
        test_spot(std::numeric_limits<T>::infinity(), fmt, i);
        test_spot(-std::numeric_limits<T>::infinity(), fmt, i);
        test_spot(std::numeric_limits<T>::quiet_NaN(), fmt, i);

        #if (defined(__clang__) && __clang_major__ >= 16 && defined(__APPLE__)) || defined(_MSC_VER)
        //
        // Newer apple clang and MSVC both give the following:
        //
        // -qNaN =  -nan(ind)
        //
        test_spot(-std::numeric_limits<T>::quiet_NaN(), fmt, i);
        #endif

        #if (defined(__clang__) && __clang_major__ >= 16 && defined(__APPLE__))
        //
        // Newer clang also gives the following:
        //
        //  sNaN =  nan(snan)
        // -sNaN = -nan(snan)
        //
        test_spot(std::numeric_limits<T>::signaling_NaN(), fmt, i);
        test_spot(-std::numeric_limits<T>::signaling_NaN(), fmt, i);
        #endif
    }
}

template <typename T>
void fixed_test()
{
    constexpr T upper_bound = std::is_same<T, double>::value ? T(std::numeric_limits<std::uint64_t>::max()) : 
                                                               T(std::numeric_limits<std::uint32_t>::max());
    
    std::mt19937_64 gen(42);
    std::uniform_real_distribution<T> dist(1, upper_bound);

    for (std::size_t i = 0; i < 1000; ++i)
    {
        for (int prec = 0; prec < 100; ++prec)
        {
            test_spot(dist(gen), boost::charconv::chars_format::fixed, prec);
        }
    }
}

// Clang 16 (most recent at time of writing) only has integral from_chars implemented
#if (defined(__GNUC__) && __GNUC__ >= 11) || (defined(_MSC_VER) && _MSC_VER >= 1924)

template <typename T>
void test_roundtrip( T value )
{
    char buffer[ 256 ];

    const auto r = std::to_chars( buffer, buffer + sizeof( buffer ), value );

    T v2 = 0;
    std::from_chars( buffer, r.ptr, v2 );

    if( BOOST_TEST_EQ( v2, value ) )
    {
    }
    else
    {
        // LCOV_EXCL_START
        #ifdef BOOST_CHARCONV_DEBUG
        std::cerr << std::setprecision(17)
                  << "     Value: " << value
                  << "\n  To chars: " << std::string( buffer, r.ptr )
                  << "\nFrom chars: " << v2 << std::endl;
        #else
        std::cerr << "... test failure for value=" << value << "; buffer='" << std::string( buffer, r.ptr ) << "'" << std::endl;
        #endif
        // LCOV_EXCL_STOP
    }
}

#endif

int main()
{   
    // General format
    random_test<float>();
    random_test<double>();
    #ifndef BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE
    random_test<long double>();
    #endif

    test_spot<double>(0.0);
    test_spot<double>(-0.0);

    // Good for debugging comparison but too many ideosyncracies
    //random_test<float>(boost::charconv::chars_format::general, -1e5F, 1e5F);
    //random_test<double>(boost::charconv::chars_format::general, -1e5, 1e5);
    //random_test<long double>(boost::charconv::chars_format::general, -1e5L, 1e5L);

    // Scientific
    random_test<float>(boost::charconv::chars_format::scientific);
    random_test<double>(boost::charconv::chars_format::scientific);
    #ifndef BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE
    random_test<long double>(boost::charconv::chars_format::scientific);
    #endif
    test_spot<double>(0.0, boost::charconv::chars_format::scientific);
    test_spot<double>(-0.0, boost::charconv::chars_format::scientific);

    random_test<float>(boost::charconv::chars_format::scientific, -1e5F, 1e5F);
    random_test<double>(boost::charconv::chars_format::scientific, -1e5, 1e5);

    // TODO(mborland): Use Ryu to fix terminal rounding rather than string manipulation
    //random_test<long double>(boost::charconv::chars_format::scientific, -1e5L, 1e5L);

    // Hex
    random_test<float>(boost::charconv::chars_format::hex);
    random_test<double>(boost::charconv::chars_format::hex);
    #ifndef BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE
    random_test<long double>(boost::charconv::chars_format::hex);
    #endif

    #if !defined(_LIBCPP_VERSION)

    random_test<float>(boost::charconv::chars_format::hex, -1e5F, 1e5F);
    random_test<double>(boost::charconv::chars_format::hex, -1e5, 1e5);
    #ifndef BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE
    random_test<long double>(boost::charconv::chars_format::hex, -1e5L, 1e5L);
    #endif

    #endif

    test_spot<double>(-9.52743282403084637e+306, boost::charconv::chars_format::hex);
    test_spot<double>(-9.52743282403084637e-306, boost::charconv::chars_format::hex);
    test_spot<double>(-9.52743282403084637e+305, boost::charconv::chars_format::hex);
    test_spot<double>(-9.52743282403084637e-305, boost::charconv::chars_format::hex);
    test_spot<double>(0.0, boost::charconv::chars_format::hex);
    test_spot<double>(-0.0, boost::charconv::chars_format::hex);

    // Fixed
    fixed_test<float>();
    fixed_test<double>();
    test_spot<double>(0.0, boost::charconv::chars_format::fixed);
    test_spot<double>(-0.0, boost::charconv::chars_format::fixed);
    test_spot<double>(0.0, boost::charconv::chars_format::fixed, 10);
    test_spot<double>(-0.0, boost::charconv::chars_format::fixed, 10);
    
    // Various non-finite values
    non_finite_test<float>();
    non_finite_test<double>();
    non_finite_test<float>(boost::charconv::chars_format::scientific);
    non_finite_test<double>(boost::charconv::chars_format::scientific);
    non_finite_test<float>(boost::charconv::chars_format::hex);
    non_finite_test<double>(boost::charconv::chars_format::hex);

    #ifndef BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE
    non_finite_test<long double>();
    non_finite_test<long double>(boost::charconv::chars_format::scientific);
    non_finite_test<long double>(boost::charconv::chars_format::hex);
    #endif

    #if (defined(__GNUC__) && __GNUC__ >= 11) || (defined(_MSC_VER) && _MSC_VER >= 1924)
    // Selected additional values
    // These are tested on boost in roundtrip.cpp
    test_roundtrip<double>(1.10393929655481808e+308);
    test_roundtrip<double>(-1.47902377240341038e+308);
    test_roundtrip<double>(-2.13177235460600904e+307);
    test_roundtrip<double>(8.60473951619578187e+307);
    test_roundtrip<double>(-2.97613696314797352e+306);

    test_roundtrip<float>(3.197633022e+38F);
    test_roundtrip<float>(2.73101834e+38F);
    test_roundtrip<float>(3.394053352e+38F);
    test_roundtrip<float>(5.549256619e+37F);
    test_roundtrip<float>(8.922125027e+34F);
    #endif

    // Reported in issue #93
    test_spot<float>(3.3F);
    test_spot<double>(3.3);

    #ifndef BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE
    // Updated tools give weird sporadic rounding errors that I can't duplicate locally
    // test_spot<long double>(3.3L);
    #endif

    return boost::report_errors();
}

#else

int main()
{
    return 0;
}

#endif
