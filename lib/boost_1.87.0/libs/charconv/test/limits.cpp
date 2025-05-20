// Copyright 2023 Matt Borland
// Copyright 2023 Peter Dimov
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

#include <boost/charconv/limits.hpp>
#include <boost/charconv/to_chars.hpp>
#include <boost/charconv/from_chars.hpp>
#include <boost/core/lightweight_test.hpp>
#include <system_error>
#include <limits>
#include <string>

void test_odr_use( int const* );

template<typename T> void test_integral( T value )
{
    // no base
    {
        char buffer[ boost::charconv::limits<T>::max_chars10 ];
        auto r = boost::charconv::to_chars( buffer, buffer + sizeof( buffer ), value );
        BOOST_TEST(r.ec == std::errc());

        T v2 = 0;
        auto r2 = boost::charconv::from_chars( buffer, r.ptr, v2 );

        if( BOOST_TEST( r2.ec == std::errc() ) && BOOST_TEST( v2 == value ) )
        {
        }
        else
        {
            std::cerr << "... test failure for value=" << value << "; buffer='" << std::string( buffer, r.ptr ) << "'" << std::endl; // LCOV_EXCL_LINE
        }
    }

    // base 10
    {
        char buffer[ boost::charconv::limits<T>::max_chars10 ];
        auto r = boost::charconv::to_chars( buffer, buffer + sizeof( buffer ), value, 10 );
        BOOST_TEST(r.ec == std::errc());

        T v2 = 0;
        auto r2 = boost::charconv::from_chars( buffer, r.ptr, v2, 10 );

        if( BOOST_TEST( r2.ec == std::errc() ) && BOOST_TEST( v2 == value ) )
        {
        }
        else
        {
            std::cerr << "... test failure for value=" << value << "; buffer='" << std::string( buffer, r.ptr ) << "'" << std::endl; // LCOV_EXCL_LINE
        }
    }

    // any base
    for( int base = 2; base <= 36; ++base )
    {
        char buffer[ boost::charconv::limits<T>::max_chars ];
        auto r = boost::charconv::to_chars( buffer, buffer + sizeof( buffer ), value, base );
        BOOST_TEST(r.ec == std::errc());

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
}

template<typename T> void test_integral()
{
    BOOST_TEST_GE( boost::charconv::limits<T>::max_chars10, std::numeric_limits<T>::digits10 );
    BOOST_TEST_GE( boost::charconv::limits<T>::max_chars, std::numeric_limits<T>::digits );

    test_odr_use( &boost::charconv::limits<T>::max_chars10 );
    test_odr_use( &boost::charconv::limits<T>::max_chars );

    test_integral( std::numeric_limits<T>::min() );
    test_integral( std::numeric_limits<T>::max() );
}

template<typename T> void test_floating_point( T value )
{
    // no base, max_chars10
    {
        char buffer[ boost::charconv::limits<T>::max_chars10 ];
        auto r = boost::charconv::to_chars( buffer, buffer + sizeof( buffer ), value );
        if (!BOOST_TEST(r.ec == std::errc()))
        {
            // LCOV_EXCL_START
            std::cerr << " Value: " << value
                      << "\nBuffer: " << std::string(buffer)
                      << "\n    Ec: " << static_cast<int>(r.ec) << std::endl;
            // LCOV_EXCL_STOP
        }

        T v2 = 0;
        auto r2 = boost::charconv::from_chars( buffer, r.ptr, v2 );

        if (!BOOST_TEST(r2.ec == std::errc()) && BOOST_TEST_EQ( v2, value ))
        {
            // LCOV_EXCL_START
            std::cerr << " Value: " << value
                      << "\nBuffer: " << std::string(buffer)
                      << "\nRetVal: " << v2
                      << "\n    Ec: " << static_cast<int>(r2.ec) << std::endl;
            // LCOV_EXCL_STOP
        }
    }

    // no base, max_chars
    {
        char buffer[ boost::charconv::limits<T>::max_chars ];
        auto r = boost::charconv::to_chars( buffer, buffer + sizeof( buffer ), value );
        if (!BOOST_TEST(r.ec == std::errc()))
        {
            // LCOV_EXCL_START
            std::cerr << " Value: " << value
                      << "\nBuffer: " << std::string(buffer)
                      << "\n    Ec: " << static_cast<int>(r.ec) << std::endl;
            // LCOV_EXCL_STOP
        }

        T v2 = 0;
        auto r2 = boost::charconv::from_chars( buffer, r.ptr, v2 );

        if (!BOOST_TEST(r2.ec == std::errc()) && BOOST_TEST_EQ( v2, value ))
        {
            // LCOV_EXCL_START
            std::cerr << " Value: " << value
                      << "\nBuffer: " << std::string(buffer)
                      << "\nRetVal: " << v2
                      << "\n    Ec: " << static_cast<int>(r2.ec) << std::endl;
            // LCOV_EXCL_STOP
        }
    }
}

template<typename T> void test_floating_point()
{
    BOOST_TEST_GE( boost::charconv::limits<T>::max_chars10, std::numeric_limits<T>::max_digits10 );
    BOOST_TEST_GE( boost::charconv::limits<T>::max_chars, std::numeric_limits<T>::max_digits10 );

    test_odr_use( &boost::charconv::limits<T>::max_chars10 );
    test_odr_use( &boost::charconv::limits<T>::max_chars );

    test_floating_point( std::numeric_limits<T>::min() );
    test_floating_point( -std::numeric_limits<T>::min() );
    test_floating_point( std::numeric_limits<T>::max() );
    test_floating_point( -std::numeric_limits<T>::max() );
}

int main()
{
    test_integral<char>();
    test_integral<signed char>();
    test_integral<unsigned char>();
    test_integral<short>();
    test_integral<unsigned short>();
    test_integral<int>();
    test_integral<unsigned>();
    test_integral<long>();
    test_integral<unsigned long>();
    test_integral<long long>();
    test_integral<unsigned long long>();

    test_floating_point<float>();
    test_floating_point<double>();

    #ifndef BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE
    test_floating_point<long double>();
    #endif

    #ifdef BOOST_CHARCONV_HAS_FLOAT16
    test_floating_point<std::float16_t>();
    #endif
    #ifdef BOOST_CHARCONV_HAS_FLOAT32
    test_floating_point<std::float32_t>();
    #endif
    #ifdef BOOST_CHARCONV_HAS_FLOAT64
    test_floating_point<std::float64_t>();
    #endif
    #ifdef BOOST_CHARCONV_HAS_BRAINFLOAT16
    test_floating_point<std::bfloat16_t>();
    #endif

#ifdef BOOST_CHARCONV_HAS_INT128

    test_integral<boost::int128_type>();
    test_integral<boost::uint128_type>();

#endif

    return boost::report_errors();
}

void test_odr_use( int const* )
{
}
