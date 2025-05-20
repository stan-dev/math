// Copyright 2022 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/config.hpp>

#ifdef BOOST_HAS_INT128

// We need to define these operator<< overloads before
// including boost/core/lightweight_test.hpp, or they
// won't be visible to BOOST_TEST_EQ
// LCOV_EXCL_START

#include <ostream>

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

// LCOV_EXCL_STOP

#endif // #ifdef BOOST_HAS_INT128

#include <boost/charconv.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/core/detail/splitmix64.hpp>
#include <system_error>
#include <iostream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <cstdint>
#include <cfloat>
#include <cmath>

int const N = 1024;

static boost::detail::splitmix64 rng;

// integral types, random values

#if defined(__GNUC__) && (__GNUC__ >= 12)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

template<class T> void test_roundtrip_integers(T value, int base )
{
    char buffer[ 256 ];

    auto r = boost::charconv::to_chars( buffer, buffer + sizeof( buffer ), value, base );

    BOOST_TEST( r.ec == std::errc() );

    T v2 = 0;
    auto r2 = boost::charconv::from_chars( buffer, r.ptr, v2, base );

    if( BOOST_TEST( r2.ec == std::errc() ) && BOOST_TEST( v2 == value ) )
    {
    }
    else
    {
        std::cerr << "... test failure for value=" << value << "; buffer='" << std::string( buffer, r.ptr ) << "'" << std::endl; // LCOV_EXCL_LINE
    }
}

#if defined(__GNUC__) && (__GNUC__ == 12)
# pragma GCC diagnostic pop
#endif

template<class T> void test_roundtrip_int8( int base )
{
    for( int i = -256; i <= 255; ++i )
    {
        test_roundtrip_integers(static_cast<T>( i ), base);
    }
}

template<class T> void test_roundtrip_uint8( int base )
{
    for( int i = 0; i <= 256; ++i )
    {
        test_roundtrip_integers(static_cast<T>( i ), base);
    }
}

template<class T> void test_roundtrip_int16( int base )
{
    test_roundtrip_int8<T>( base );

    for( int i = 0; i < N; ++i )
    {
        std::int16_t w = static_cast<std::int16_t>( rng() );
        test_roundtrip_integers(static_cast<T>( w ), base);
    }
}

template<class T> void test_roundtrip_uint16( int base )
{
    test_roundtrip_uint8<T>( base );

    for( int i = 0; i < N; ++i )
    {
        std::uint16_t w = static_cast<std::uint16_t>( rng() );
        test_roundtrip_integers(static_cast<T>( w ), base);
    }
}

template<class T> void test_roundtrip_int32( int base )
{
    test_roundtrip_int16<T>( base );

    for( int i = 0; i < N; ++i )
    {
        std::int32_t w = static_cast<std::int32_t>( rng() );
        test_roundtrip_integers(static_cast<T>( w ), base);
    }
}

template<class T> void test_roundtrip_uint32( int base )
{
    test_roundtrip_uint16<T>( base );

    for( int i = 0; i < N; ++i )
    {
        std::uint32_t w = static_cast<std::uint32_t>( rng() );
        test_roundtrip_integers(static_cast<T>( w ), base);
    }
}

template<class T> void test_roundtrip_int64( int base )
{
    test_roundtrip_int32<T>( base );

    for( int i = 0; i < N; ++i )
    {
        std::int64_t w = static_cast<std::int64_t>( rng() );
        test_roundtrip_integers(static_cast<T>( w ), base);
    }
}

template<class T> void test_roundtrip_uint64( int base )
{
    test_roundtrip_uint32<T>( base );

    for( int i = 0; i < N; ++i )
    {
        std::uint64_t w = static_cast<std::uint64_t>( rng() );
        test_roundtrip_integers(static_cast<T>( w ), base);
    }
}

#ifdef BOOST_CHARCONV_HAS_INT128

inline boost::uint128_type concatenate(std::uint64_t word1, std::uint64_t word2)
{
    return static_cast<boost::uint128_type>(word1) << 64 | word2;
}

template<class T> void test_roundtrip_int128( int base )
{
    for( int i = 0; i < N; ++i )
    {
        boost::int128_type w = static_cast<boost::int128_type>( concatenate(rng(), rng()) );
        test_roundtrip_integers(static_cast<T>( w ), base);
    }
}

template<class T> void test_roundtrip_uint128( int base )
{
    for( int i = 0; i < N; ++i )
    {
        boost::uint128_type w = static_cast<boost::uint128_type>( concatenate(rng(), rng()) );
        test_roundtrip_integers(static_cast<T>( w ), base);
    }
}

#endif // #ifdef BOOST_CHARCONV_HAS_INT128

// integral types, boundary values

template<class T> void test_roundtrip_bv( int base )
{
    test_roundtrip_integers(std::numeric_limits<T>::min(), base);
    test_roundtrip_integers(std::numeric_limits<T>::max(), base);
}

#ifdef BOOST_CHARCONV_HAS_INT128
template <> void test_roundtrip_bv<boost::int128_type>(int base)
{
    test_roundtrip_integers(BOOST_CHARCONV_INT128_MIN, base);
    test_roundtrip_integers(BOOST_CHARCONV_INT128_MAX, base);
}

template <> void test_roundtrip_bv<boost::uint128_type>(int base)
{
    test_roundtrip_integers(0, base);
    test_roundtrip_integers(BOOST_CHARCONV_UINT128_MAX, base);
}
#endif

// floating point types, random values

template<class T> void test_roundtrip( T value, boost::charconv::chars_format fmt = boost::charconv::chars_format::general )
{
    char buffer[ 256 ];

    auto r = boost::charconv::to_chars( buffer, buffer + sizeof( buffer ), value, fmt );

    BOOST_TEST( r.ec == std::errc() );

    T v2 = 0;
    auto r2 = boost::charconv::from_chars( buffer, r.ptr, v2, fmt );

    if( BOOST_TEST( r2.ec == std::errc() ) && BOOST_TEST_EQ( v2, value ) && BOOST_TEST( r2.ptr == r.ptr) )
    {
    }
    else
    {
        // LCOV_EXCL_START
        #ifdef BOOST_CHARCONV_DEBUG_ROUNDTRIP
        std::cerr << std::setprecision(std::numeric_limits<T>::digits10 + 1)
                  << "     Value: " << value
                  << "\n  To chars: " << std::string( buffer, r.ptr )
                  << "\nFrom chars: " << v2 << std::endl
                  << std::hexfloat
                  << "\n     Value: " << value
                  << "\nFrom chars: " << v2
                  << "\n R1 offset: " << (r.ptr - buffer)
                  << "\n R2 offset: " << (r2.ptr - buffer) << std::endl << std::scientific;
        #else
        std::cerr << "... test failure for value=" << value << "; buffer='" << std::string( buffer, r.ptr ) << "'" << std::endl;
        #endif
        // LCOV_EXCL_STOP
    }
}

// https://stackoverflow.com/questions/62074229/float-distance-for-80-bit-long-double
/*  Return the signed distance from 0 to x, measuring distance as one unit per
    number representable in FPType.  x must be a finite number.
*/

#if defined(__GNUC__) && (__GNUC__ >= 5)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wattributes"
#endif

template<typename FPType>
#if !defined(BOOST_MSVC) && !(defined(__clang__) && (__clang_major__ == 3) && (__clang_minor__ < 7))
__attribute__((no_sanitize("undefined")))
#endif
int64_t ToOrdinal(FPType x)
{
    static constexpr int
            Radix             = std::numeric_limits<FPType>::radix,
            SignificandDigits = std::numeric_limits<FPType>::digits,
            MinimumExponent   = std::numeric_limits<FPType>::min_exponent;

    //  Number of normal representable numbers for each exponent.
    static const auto
            NumbersPerExponent = static_cast<uint64_t>(scalbn(Radix-1, SignificandDigits-1));

    if (x == 0)
        return 0;

    //  Separate the sign.
    int sign = std::signbit(x) ? -1 : +1;
    x = std::fabs(x);

    //  Separate the significand and exponent.
    int exponent = std::ilogb(x)+1;
    FPType fraction = std::scalbn(x, -exponent);

    if (exponent < MinimumExponent)
    {
        //  For subnormal x, adjust to its subnormal representation.
        fraction = std::scalbn(fraction, exponent - MinimumExponent);
        exponent = MinimumExponent;
    }

    /*  Start with the number of representable numbers in preceding normal
        exponent ranges.
    */
    auto count = static_cast<int64_t>(static_cast<uint64_t>(exponent - MinimumExponent) * NumbersPerExponent);

    /*  For subnormal numbers, fraction * radix ** SignificandDigits is the
        number of representable numbers from 0 to x.  For normal numbers,
        (fraction-1) * radix ** SignificandDigits is the number of
        representable numbers from the start of x's exponent range to x, and
        1 * radix ** SignificandDigits is the number of representable subnormal
        numbers (which we have not added into count yet).  So, in either case,
        adding fraction * radix ** SignificandDigits is the desired amount to
        add to count.
    */
    count += static_cast<int64_t>(std::scalbn(fraction, SignificandDigits));

    return sign * count;
}

#if defined(__GNUC__) && (__GNUC__ >= 5)
# pragma GCC diagnostic pop
#endif

/*  Return the number of representable numbers from x to y, including one
    endpoint.
*/
template<typename FPType> int64_t Distance(FPType y, FPType x)
{
    return ToOrdinal(y) - ToOrdinal(x);
}

#ifndef BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE
template <> void test_roundtrip<long double>(long double value, boost::charconv::chars_format fmt)
{
    char buffer[ 256 ];

    auto r = boost::charconv::to_chars( buffer, buffer + sizeof( buffer ), value, fmt );

    BOOST_TEST( r.ec == std::errc() );

    long double v2 = 0;
    auto r2 = boost::charconv::from_chars( buffer, r.ptr, v2, fmt );

    if( BOOST_TEST( r2.ec == std::errc() ) && BOOST_TEST( std::abs(Distance(v2, value)) < INT64_C(1) ) )
    {
    }
    else
    {
        // LCOV_EXCL_START
        #ifdef BOOST_CHARCONV_DEBUG_ROUNDTRIP
        std::cerr << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
                  << "     Value: " << value
                  << "\n  To chars: " << std::string( buffer, r.ptr )
                  << "\nFrom chars: " << v2 << std::endl
                  << std::hexfloat
                  << "\n     Value: " << value
                  << "\nFrom chars: " << v2 << std::endl << std::scientific;
        #else
        std::cerr << "... test failure for value=" << value
                  << "; buffer='" << std::string( buffer, r.ptr ) << "'"
                  << "; ulp distance=" << Distance(v2, value)
                  << "; error code=" << static_cast<int>(r2.ec) << std::endl;
        #endif
        // LCOV_EXCL_STOP
    }
}
#endif

// floating point types, boundary values

template<class T> void test_roundtrip_bv()
{
    test_roundtrip( std::numeric_limits<T>::min() );
    test_roundtrip( -std::numeric_limits<T>::min() );
    test_roundtrip( std::numeric_limits<T>::max() );
    test_roundtrip( +std::numeric_limits<T>::max() );
}

//

template <typename T>
void test_extreme_values()
{
    T current_pos = (std::numeric_limits<T>::max)();
    for (int i = 0; i < 10000; ++i)
    {
        test_roundtrip<T>(current_pos);
        current_pos = std::nexttoward(current_pos, T(0));
    }

    current_pos = (std::numeric_limits<T>::min)();
    for (int i = 0; i < 10000; ++i)
    {
        test_roundtrip<T>(current_pos);
        current_pos = std::nextafter(current_pos, T(1));
    }
}

int main()
{
    // integral types, random values

    for( int base = 2; base <= 36; ++base )
    {
        test_roundtrip_int8<std::int8_t>( base );
        test_roundtrip_uint8<std::uint8_t>( base );

        test_roundtrip_int16<std::int16_t>( base );
        test_roundtrip_uint16<std::uint16_t>( base );

        test_roundtrip_int32<std::int32_t>( base );
        test_roundtrip_uint32<std::uint32_t>( base );

        test_roundtrip_int64<std::int64_t>( base );
        test_roundtrip_uint64<std::uint64_t>( base );

#ifdef BOOST_CHARCONV_HAS_INT128

        test_roundtrip_int128<boost::int128_type>( base );
        test_roundtrip_uint128<boost::uint128_type>( base );

#endif
    }

    // integral types, boundary values

    for( int base = 2; base <= 36; ++base )
    {
        test_roundtrip_bv<char>( base );
        test_roundtrip_bv<signed char>( base );
        test_roundtrip_bv<unsigned char>( base );

        test_roundtrip_bv<short>( base );
        test_roundtrip_bv<unsigned short>( base );

        test_roundtrip_bv<int>( base );
        test_roundtrip_bv<unsigned int>( base );

        test_roundtrip_bv<long>( base );
        test_roundtrip_bv<unsigned long>( base );

        test_roundtrip_bv<long long>( base );
        test_roundtrip_bv<unsigned long long>( base );

#ifdef BOOST_CHARCONV_HAS_INT128

        test_roundtrip_bv<boost::int128_type>( base );
        test_roundtrip_bv<boost::uint128_type>( base );

#endif
    }

    // 16-bit types

    double const q = std::pow( 1.0, -64 );

    #ifdef BOOST_CHARCONV_HAS_FLOAT16
    {
        std::float16_t const small_q = std::pow(1.0F16, -16.0F16);

        for( int i = 0; i < N; ++i )
        {
            std::float16_t w0 = static_cast<std::float16_t>( rng() ); // 0 .. 2^64
            test_roundtrip( w0 );

            std::float16_t w1 = static_cast<std::float16_t>( rng() ) * small_q ; // 0.0 .. 1.0
            test_roundtrip( w1 );

            std::float16_t w2 = (std::numeric_limits<std::float16_t>::max)() / static_cast<std::float16_t>( rng() ); // large values
            test_roundtrip( w2 );

            std::float16_t w3 = (std::numeric_limits<std::float16_t>::min)() * static_cast<std::float16_t>( rng() ); // small values
            test_roundtrip( w3 );
        }

        test_roundtrip_bv<std::float16_t>();
    }
    #endif

    #ifdef BOOST_CHARCONV_HAS_BRAINFLOAT16
    {
        std::bfloat16_t const small_q = std::pow(1.0BF16, -16.0BF16);

        for( int i = 0; i < N; ++i )
        {
            std::bfloat16_t w0 = static_cast<std::bfloat16_t>( rng() ); // 0 .. 2^64
            test_roundtrip( w0 );

            std::bfloat16_t w1 = static_cast<std::bfloat16_t>( rng() ) * small_q ; // 0.0 .. 1.0
            test_roundtrip( w1 );

            std::bfloat16_t w2 = (std::numeric_limits<std::bfloat16_t>::max)() / static_cast<std::bfloat16_t>( rng() ); // large values
            test_roundtrip( w2 );

            std::bfloat16_t w3 = (std::numeric_limits<std::bfloat16_t>::min)() * static_cast<std::bfloat16_t>( rng() ); // small values
            test_roundtrip( w3 );
        }

        test_roundtrip_bv<std::bfloat16_t>();
    }
    #endif

    // float

    {
        for( int i = 0; i < N; ++i )
        {
            float w0 = static_cast<float>( rng() ); // 0 .. 2^64
            test_roundtrip( w0 );
            test_roundtrip( w0, boost::charconv::chars_format::fixed );

            float w1 = static_cast<float>( rng() ) * static_cast<float>( q ); // 0.0 .. 1.0
            test_roundtrip( w1 );
            test_roundtrip( w1, boost::charconv::chars_format::fixed );

            float w2 = FLT_MAX / static_cast<float>( rng() ); // large values
            test_roundtrip( w2 );
            test_roundtrip( w2, boost::charconv::chars_format::fixed );

            float w3 = FLT_MIN * static_cast<float>( rng() ); // small values
            test_roundtrip( w3 );
            test_roundtrip( w3, boost::charconv::chars_format::fixed );
        }

        test_roundtrip_bv<float>();
    }

    #ifdef BOOST_CHARCONV_HAS_FLOAT32
    {
        for( int i = 0; i < N; ++i )
        {
            std::float32_t w0 = static_cast<std::float32_t>( rng() ); // 0 .. 2^64
            test_roundtrip( w0 );

            std::float32_t w1 = static_cast<std::float32_t>( rng() ) *  static_cast<std::float32_t>(q) ; // 0.0 .. 1.0
            test_roundtrip( w1 );

            std::float32_t w2 = (std::numeric_limits<std::float32_t>::max)() / static_cast<std::float32_t>( rng() ); // large values
            test_roundtrip( w2 );

            std::float32_t w3 = (std::numeric_limits<std::float32_t>::min)() * static_cast<std::float32_t>( rng() ); // small values
            test_roundtrip( w3 );
        }

        test_roundtrip_bv<std::float32_t>();
    }
    #endif

    // double

    {
        for( int i = 0; i < N; ++i )
        {
            double w0 = static_cast<double>( rng() ) * 1.0; // 0 .. 2^64
            test_roundtrip( w0 );
            test_roundtrip( w0, boost::charconv::chars_format::fixed );

            double w1 = static_cast<double>( rng() ) * q; // 0.0 .. 1.0
            test_roundtrip( w1 );
            test_roundtrip( w1, boost::charconv::chars_format::fixed );

            double w2 = DBL_MAX / static_cast<double>( rng() ); // large values
            test_roundtrip( w2 );

            double w3 = DBL_MIN * static_cast<double>( rng() ); // small values
            test_roundtrip( w3 );
        }

        test_roundtrip_bv<double>();
    }

    #ifdef BOOST_CHARCONV_HAS_FLOAT64
    {
        for( int i = 0; i < N; ++i )
        {
            std::float64_t w0 = static_cast<std::float64_t>( rng() ); // 0 .. 2^64
            test_roundtrip( w0 );

            std::float64_t w1 = static_cast<std::float64_t>( rng() ) * static_cast<std::float64_t>(q) ; // 0.0 .. 1.0
            test_roundtrip( w1 );

            std::float64_t w2 = (std::numeric_limits<std::float64_t>::max)() / static_cast<std::float64_t>( rng() ); // large values
            test_roundtrip( w2 );

            std::float64_t w3 = (std::numeric_limits<std::float64_t>::min)() * static_cast<std::float64_t>( rng() ); // small values
            test_roundtrip( w3 );
        }

        test_roundtrip_bv<std::float64_t>();
    }
    #endif

    // long double
    #if !(BOOST_CHARCONV_LDBL_BITS == 128) && !defined(BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE)

    {
        long double const ql = std::pow( 1.0L, -64 );

        for( int i = 0; i < N; ++i )
        {
            long double w0 = static_cast<long double>( rng() ) * 1.0L; // 0 .. 2^64
            test_roundtrip( w0 );

            long double w1 = static_cast<long double>( rng() ) * ql; // 0.0 .. 1.0
            test_roundtrip( w1 );

            long double w2 = LDBL_MAX / static_cast<long double>( rng() ); // large values
            test_roundtrip( w2 );

            long double w3 = LDBL_MIN * static_cast<long double>( rng() ); // small values
            test_roundtrip( w3 );

            long double w4 = -static_cast<long double>( rng() ) * 1.0L; // -0 .. 2^64
            test_roundtrip( w4 );
        }

        test_roundtrip_bv<long double>();
    }

    #endif

    // Selected additional values
    //
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

    test_extreme_values<float>();
    test_extreme_values<double>();

    return boost::report_errors();
}
