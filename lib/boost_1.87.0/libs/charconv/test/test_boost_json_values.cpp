// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
// Copyright (c) 2020 Krystian Stasiowski (sdkrystian@gmail.com)
// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/json/issues/599
// See: https://github.com/boostorg/json/blob/develop/test/double.cpp

#include <boost/charconv.hpp>
#include <boost/core/lightweight_test.hpp>
#include <system_error>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <limits>
#include <random>
#include <type_traits>
#include <cstring>
#include <cinttypes>

template <typename T>
void grind(const std::string& str, const T expected_value)
{
    // From string to expected value
    T v {};
    auto from_r = boost::charconv::from_chars(str.c_str(), str.c_str() + std::strlen(str.c_str()), v);
    if (!(BOOST_TEST_EQ(v, expected_value) && BOOST_TEST(from_r.ec == std::errc())))
    {
        // LCOV_EXCL_START
        std::cerr << "Expected value: " << expected_value << "\nFrom chars value: " << v << std::endl;
        return;
        // LCOV_EXCL_STOP
    }

    // Roundtrip
    T roundtrip_v {};
    char buffer[256] {};
    auto to_r = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), expected_value);
    BOOST_TEST(to_r.ec == std::errc());

    auto roundtrip_r = boost::charconv::from_chars(buffer, buffer + std::strlen(buffer), roundtrip_v);
    if (!(BOOST_TEST_EQ(roundtrip_v, expected_value) && BOOST_TEST(roundtrip_r.ec == std::errc())))
    {
        // LCOV_EXCL_START
        std::cerr << "Expected value: " << expected_value << "\nRoundtrip value: " << roundtrip_v << std::endl;
        return;
        // LCOV_EXCL_STOP
    }
}

void grind_double(const std::string& str, const double expected_value)
{
    grind<double>(str, expected_value);
}

void grind_int64(const std::string& str, const std::int64_t expected_value)
{
    grind<std::int64_t>(str, expected_value);
}

void grind_uint64(const std::string& str, const std::uint64_t expected_value)
{
    grind<std::uint64_t>(str, expected_value);
}

template <typename T>
void spot_value(const std::string& buffer, T expected_value)
{
    T v = 0;
    auto r = boost::charconv::from_chars_erange(buffer.c_str(), buffer.c_str() + std::strlen(buffer.c_str()), v);
    BOOST_TEST(r.ec == std::errc());
    if (!BOOST_TEST_EQ(v, expected_value))
    {
        std::cerr << "Test failure for: " << buffer << " got: " << v << std::endl; // LCOV_EXCL_LINE
    }
}

void fc (const std::string& s)
{
    char* str_end;
    const double expected_value = std::strtod(s.c_str(), &str_end);
    spot_value(s, expected_value);
}

// https://github.com/boostorg/json/issues/599
void issue_599_test()
{
    const std::vector<double> ref_values = {
            -0.27006296145688152
            , -0.42448451824686562
            , 0.057166336253381224
            , 1217.2772861138403
            , -161.68713249779881
            , 267.04251495962637
            , -0.66615716744854903
            , -0.80918535242172396
            , 0.29123256908034584
            , 2137.0241926849581
            , -599.61476423470071
            , 1002.9111801605201
            , -1.4239725866123634
            , -1.0346132515963697
            , 0.90790798114618365
            , 2535.2404973980229
            , -1207.1290929173115
            , 2028.379668845469
            , -2.2996584528580817
            , -0.90521880279912548
            , 1.6449616573025234
            , 2314.9160231072947
            , -1836.2705293282695
            , 3093.6759124836558
            , -2.7781953227473752
            , -0.54944860097807424
            , 1.9702410871568496
            , 1529.7366713650281
            , -2333.8061352785221
            , 3939.3529821428001
            , -3.0696620243882053
            , -0.13139419714747352
            , 2.0689802496328444
            , 370.91535570754445
            , -2578.5523049665962
            , 4359.4347976947429
            , 2.9538186832340876
            , 0.29606564516163103
            , 2.0456322162820761
            , -879.28776788268817
            , -2510.8913760620435
            , 4251.6098535812462 };

    for (const auto current_ref_val : ref_values)
    {
        char buffer[256] {};
        const auto r = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), current_ref_val);
        BOOST_TEST(r.ec == std::errc());

        double return_val {};
        const auto return_r = boost::charconv::from_chars_erange(buffer, buffer + std::strlen(buffer), return_val);
        BOOST_TEST(return_r.ec == std::errc());
        if (!BOOST_TEST_EQ(current_ref_val, return_val))
        {
            // LCOV_EXCL_START
            #ifdef BOOST_CHARCONV_DEBUG
            std::cerr << std::setprecision(17)
                  << "     Value: " << current_ref_val
                  << "\n  To chars: " << std::string( buffer, r.ptr )
                  << "\nFrom chars: " << return_val << std::endl;
            #else
            std::cerr << "... test failure for value=" << current_ref_val << "; buffer='" << std::string( buffer, r.ptr ) << "'" << std::endl;
            #endif
            // LCOV_EXCL_STOP
        }
    }
}

#ifdef BOOST_MSVC
# pragma warning(push)
# pragma warning(disable: 4244) // Conversion from double to T with BOOST_IF_CONSTEXPR expansion pre-C++17
#elif defined(__GNUC__) && __GNUC__ >= 5
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wmissing-field-initializers"
# pragma GCC diagnostic ignored "-Wfloat-conversion"
#elif defined(__clang__)
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wfloat-conversion"
# pragma clang diagnostic ignored "-Wconversion"
# if __clang_major__ > 7
#   pragma clang diagnostic ignored "-Wimplicit-float-conversion"
# endif
#endif

template <typename T>
void check_accuracy(const char* nm, int max_ulp)
{
    using Unsigned_Integer = typename std::conditional<std::is_same<T, float>::value, std::uint32_t, std::uint64_t>::type;

    T x {};
    BOOST_IF_CONSTEXPR (std::is_same<T, float>::value)
    {
        x = std::strtof( nm, nullptr );
    }
    else
    {
        x = std::strtod( nm, nullptr );
    }

    T y {};
    boost::charconv::from_chars_erange(nm, nm + std::strlen(nm), y, boost::charconv::chars_format::scientific);

    Unsigned_Integer bx;
    Unsigned_Integer by;
    std::memcpy( &bx, &x, sizeof(x) );
    std::memcpy( &by, &y, sizeof(y) );
    const auto diff = static_cast<std::int64_t>(bx - by);
    if (!BOOST_TEST(std::abs( diff ) <= max_ulp))
        std::fprintf(stderr,
                     "%s: difference %" PRId64 " ulp\n"
                                               "  strtod:       %.13a %.16g\n"
                                               "  boost.json:   %.13a %.16g\n\n",
            nm, diff, static_cast<double>(x), static_cast<double>(x), static_cast<double>(y), static_cast<double>(y) );
}

#ifdef BOOST_MSVC
# pragma warning(pop)
#elif defined(__GNUC__) && __GNUC__ >= 5
# pragma GCC diagnostic pop
#elif defined(__clang__)
# pragma clang diagnostic pop
#endif

template <typename T>
void test_within_ulp()
{
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<> dist( std::numeric_limits<T>::min_exponent10 - 1, std::numeric_limits<T>::max_exponent10 + 1);

    check_accuracy<T>("10199214983525025199.13135016100190689227e-308", 2);

    for( int i = 0; i < 1000000; ++i )
    {
        unsigned long long x1 = rng();
        unsigned long long x2 = rng();
        int x3 = dist( rng );

        char buffer[ 128 ];
        snprintf( buffer, sizeof(buffer), "%llu.%llue%d", x1, x2, x3 );

        check_accuracy<T>(buffer, 1);
    }

    for( int i = -326; i <= +309; ++i )
    {
        char buffer[ 128 ];
        snprintf( buffer, sizeof(buffer), "1e%d", i );

        check_accuracy<T>(buffer, 0);
    }
};

int main()
{
    issue_599_test();

    grind_double("-1.010", -1.01);
    grind_double("-0.010", -0.01);
    grind_double("-0.0", -0.0);
    grind_double("-0e0", -0.0);
    grind_double("18.4",  18.4);
    grind_double("-18.4", -18.4);
    grind_double("18446744073709551616",  1.8446744073709552e+19);
    grind_double("-18446744073709551616", -1.8446744073709552e+19);
    grind_double("18446744073709551616.0",  1.8446744073709552e+19);
    grind_double("18446744073709551616.00009",  1.8446744073709552e+19);
    grind_double("1844674407370955161600000",  1.8446744073709552e+24);
    grind_double("-1844674407370955161600000", -1.8446744073709552e+24);
    grind_double("1844674407370955161600000.0",  1.8446744073709552e+24);
    grind_double("1844674407370955161600000.00009",  1.8446744073709552e+24);
    grind_double("19700720435664.186294290058937593e13",  1.9700720435664185e+26);

    grind_double("1.0", 1.0);
    grind_double("1.1", 1.1);
    grind_double("1.11", 1.11);
    grind_double("1.11111", 1.11111);
    grind_double("11.1111", 11.1111);
    grind_double("111.111", 111.111);

    // From https://github.com/boostorg/json/blob/develop/test/stream_parser.cpp#L673
    grind_int64( "-9223372036854775808", INT64_MIN);
    grind_int64( "-9223372036854775807", -9223372036854775807);
    grind_int64( "-999999999999999999", -999999999999999999);
    grind_int64( "-99999999999999999", -99999999999999999);
    grind_int64( "-9999999999999999", -9999999999999999);
    grind_int64( "-999999999999999", -999999999999999);
    grind_int64( "-99999999999999", -99999999999999);
    grind_int64( "-9999999999999", -9999999999999);
    grind_int64( "-999999999999", -999999999999);
    grind_int64( "-99999999999", -99999999999);
    grind_int64( "-9999999999", -9999999999);
    grind_int64( "-999999999", -999999999);
    grind_int64( "-99999999", -99999999);
    grind_int64( "-9999999", -9999999);
    grind_int64( "-999999", -999999);
    grind_int64( "-99999", -99999);
    grind_int64( "-9999", -9999);
    grind_int64( "-999", -999);
    grind_int64( "-99", -99);
    grind_int64( "-9", -9);
    grind_int64( "-0", 0);
    grind_int64(  "0", 0);
    grind_int64(  "1", 1);
    grind_int64(  "9", 9);
    grind_int64(  "99", 99);
    grind_int64(  "999", 999);
    grind_int64(  "9999", 9999);
    grind_int64(  "99999", 99999);
    grind_int64(  "999999", 999999);
    grind_int64(  "9999999", 9999999);
    grind_int64(  "99999999", 99999999);
    grind_int64(  "999999999", 999999999);
    grind_int64(  "9999999999", 9999999999);
    grind_int64(  "99999999999", 99999999999);
    grind_int64(  "999999999999", 999999999999);
    grind_int64(  "9999999999999", 9999999999999);
    grind_int64(  "99999999999999", 99999999999999);
    grind_int64(  "999999999999999", 999999999999999);
    grind_int64(  "9999999999999999", 9999999999999999);
    grind_int64(  "99999999999999999", 99999999999999999);
    grind_int64(  "999999999999999999", 999999999999999999);
    grind_int64(  "9223372036854775807", INT64_MAX);

    grind_uint64( "9223372036854775808", 9223372036854775808ULL);
    grind_uint64( "9999999999999999999", 9999999999999999999ULL);
    grind_uint64( "18446744073709551615", UINT64_MAX);

    grind_double("-1.010", -1.01);
    grind_double("-0.010", -0.01);
    grind_double("-0.0", -0.0);
    grind_double("-0e0", -0.0);
    grind_double( "18446744073709551616",  1.8446744073709552e+19);
    grind_double("-18446744073709551616", -1.8446744073709552e+19);
    grind_double( "18446744073709551616.0",  1.8446744073709552e+19);
    grind_double( "18446744073709551616.00009",  1.8446744073709552e+19);
    grind_double( "1844674407370955161600000",  1.8446744073709552e+24);
    grind_double("-1844674407370955161600000", -1.8446744073709552e+24);
    grind_double( "1844674407370955161600000.0",  1.8446744073709552e+24);
    grind_double( "1844674407370955161600000.00009",  1.8446744073709552e+24);

    grind_double( "1.0", 1.0);
    grind_double( "1.1", 1.1);
    grind_double( "1.11", 1.11);
    grind_double( "1.11111", 1.11111);
    grind_double( "11.1111", 11.1111);
    grind_double( "111.111", 111.111);

    fc("-999999999999999999999");
    fc("-100000000000000000009");
    fc("-10000000000000000000");
    fc("-9223372036854775809");

    fc("18446744073709551616");
    fc("99999999999999999999");
    fc("999999999999999999999");
    fc("1000000000000000000000");
    fc("9999999999999999999999");
    fc("99999999999999999999999");

    fc("-0.9999999999999999999999");
    fc("-0.9999999999999999");
    fc("-0.9007199254740991");
    fc("-0.999999999999999");
    fc("-0.99999999999999");
    fc("-0.9999999999999");
    fc("-0.999999999999");
    fc("-0.99999999999");
    fc("-0.9999999999");
    fc("-0.999999999");
    fc("-0.99999999");
    fc("-0.9999999");
    fc("-0.999999");
    fc("-0.99999");
    fc("-0.9999");
    fc("-0.8125");
    fc("-0.999");
    fc("-0.99");
    fc("-1.0");
    fc("-0.9");
    fc("-0.0");
    fc("0.0");
    fc("0.9");
    fc("0.99");
    fc("0.999");
    fc("0.8125");
    fc("0.9999");
    fc("0.99999");
    fc("0.999999");
    fc("0.9999999");
    fc("0.99999999");
    fc("0.999999999");
    fc("0.9999999999");
    fc("0.99999999999");
    fc("0.999999999999");
    fc("0.9999999999999");
    fc("0.99999999999999");
    fc("0.999999999999999");
    fc("0.9007199254740991");
    fc("0.9999999999999999");
    fc("0.9999999999999999999999");
    fc("0.999999999999999999999999999");

    fc("-1e308");
    fc("-1e-308");
    fc("-9999e300");
    fc("-999e100");
    fc("-99e10");
    fc("-9e1");
    fc("9e1");
    fc("99e10");
    fc("999e100");
    fc("9999e300");
    fc("999999999999999999.0");
    fc("999999999999999999999.0");
    fc("999999999999999999999e5");
    fc("999999999999999999999.0e5");

    fc("0.00000000000000001");

    fc("-1e-1");
    fc("-1e0");
    fc("-1e1");
    fc("0e0");
    fc("1e0");
    fc("1e10");

    fc("0."
       "00000000000000000000000000000000000000000000000000" // 50 zeroes
       "1e50");
    fc("-0."
       "00000000000000000000000000000000000000000000000000" // 50 zeroes
       "1e50");

    fc("0."
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000" // 500 zeroes
       "1e600");

    fc("-0."
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000" // 500 zeroes
       "1e600");

    fc("0e"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000"
       "00000000000000000000000000000000000000000000000000" // 500 zeroes
    );

    // Reported in issue #29
    // https://github.com/cppalliance/charconv/issues/29

    fc("0E0");
    fc("0E01");
    fc("0.0e0");
    fc("-0E0");
    fc("-0E01");
    fc("-0.0e0");

    // Reported in issue #39
    // https://github.com/cppalliance/charconv/issues/39
    fc("13037152512515243620.9123444678836069176e0");
    fc("17697431460539906388.13521888826860801885e0");
    fc("17947220026488055965.10643827068131766968e0");
    fc("11270160277506678000.2887804927902173715e-2");
    fc("16014537415383412664.8126122414435723949e0");
    test_within_ulp<double>();
    test_within_ulp<float>();
    #ifdef BOOST_CHARCONV_FULL_LONG_DOUBLE_IMPL
    test_within_ulp<long double>();
    #endif

    return boost::report_errors();
}
