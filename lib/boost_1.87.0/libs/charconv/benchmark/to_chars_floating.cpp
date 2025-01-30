// Copyright 2023 Peter Dimov
// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/charconv/detail/config.hpp>
#include <ostream>

#ifdef BOOST_CHARCONV_HAS_STDFLOAT128
#include <charconv>

std::ostream& operator<<( std::ostream& os, std::float128_t v)
{
    char buffer [ 256 ] {};
    std::to_chars(buffer, buffer + sizeof(buffer), v);
    os << buffer;
    return os;
}

#endif


#include <boost/charconv/to_chars.hpp>
#include <boost/core/type_name.hpp>
#include <boost/core/detail/splitmix64.hpp>
#include <boost/config.hpp>
#include <type_traits>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <charconv>

constexpr unsigned N = 2'000'000;
constexpr int K = 10;

template<class T> static BOOST_NOINLINE void init_input_data( std::vector<T>& data )
{
    data.reserve( N );

    boost::detail::splitmix64 rng;

    for( unsigned i = 0; i < N; ++i )
    {
        std::uint64_t tmp = rng();

        T x;
        std::memcpy( &x, &tmp, sizeof(x) );

        if( !std::isfinite(x) ) continue;

        data.push_back( x );
    }
}

#ifdef BOOST_CHARCONV_HAS_STDFLOAT128
template<> BOOST_NOINLINE void init_input_data<std::float128_t>( std::vector<std::float128_t>& data )
{
    data.reserve( N );

    boost::detail::splitmix64 rng;

    for( unsigned i = 0; i < N; ++i )
    {
        boost::charconv::detail::uint128 tmp {rng(), rng()};
        boost::uint128_type temp {tmp};

        std::float128_t x;
        std::memcpy( &x, &temp, sizeof(x) );

        if( !std::isfinite(x) ) continue;

        data.push_back( x );
    }
}
#endif

using namespace std::chrono_literals;

template<class T> static BOOST_NOINLINE void test_snprintf( std::vector<T> const& data, bool general, char const* label, int precision )
{
    auto t1 = std::chrono::steady_clock::now();

    std::size_t s = 0;

    char const* format = general? "%.*g": "%.*e";

    int prec = precision;

    if( prec == 0 )
    {
        prec = 16;

        BOOST_IF_CONSTEXPR( std::is_same<T, float>::value )
        {
            prec = 8;
        }
    }

    for( int i = 0; i < K; ++i )
    {
        char buffer[ 22 ];

        for( auto x: data )
        {
            auto r = std::snprintf( buffer, sizeof( buffer ), format, prec, x );
            s += r;
            s += static_cast<unsigned char>( buffer[0] );
        }
    }

    auto t2 = std::chrono::steady_clock::now();

    std::cout << "            std::snprintf<" << boost::core::type_name<T>() << ">, " << label << ", " << precision << ": " << std::setw( 5 ) << ( t2 - t1 ) / 1ms << " ms (s=" << s << ")\n";
}

template<class T> static BOOST_NOINLINE void test_std_to_chars( std::vector<T> const& data, bool general, char const* label, int precision )
{
    auto t1 = std::chrono::steady_clock::now();

    std::size_t s = 0;

    for( int i = 0; i < K; ++i )
    {
        char buffer[ 21 ];

        for( auto x: data )
        {
            std::chars_format fmt = general? std::chars_format::general: std::chars_format::scientific;

            auto r = precision == 0?
                     std::to_chars( buffer, buffer + sizeof( buffer ), x, fmt ):
                     std::to_chars( buffer, buffer + sizeof( buffer ), x, fmt, precision );

            s += static_cast<std::size_t>( r.ptr - buffer );
            s += static_cast<unsigned char>( buffer[0] );
        }
    }

    auto t2 = std::chrono::steady_clock::now();

    std::cout << "            std::to_chars<" << boost::core::type_name<T>() << ">, " << label << ", " << precision << ": " << std::setw( 5 ) << ( t2 - t1 ) / 1ms << " ms (s=" << s << ")\n";
}

template<class T> static BOOST_NOINLINE void test_boost_to_chars( std::vector<T> const& data, bool general, char const* label, int precision )
{
    auto t1 = std::chrono::steady_clock::now();

    std::size_t s = 0;

    for( int i = 0; i < K; ++i )
    {
        char buffer[ 21 ];

        for( auto x: data )
        {
            boost::charconv::chars_format fmt = general? boost::charconv::chars_format::general: boost::charconv::chars_format::scientific;

            auto r = precision == 0?
                     boost::charconv::to_chars( buffer, buffer + sizeof( buffer ), x, fmt ):
                     boost::charconv::to_chars( buffer, buffer + sizeof( buffer ), x, fmt, precision );

            s += static_cast<std::size_t>( r.ptr - buffer );
            s += static_cast<unsigned char>( buffer[0] );
        }
    }

    auto t2 = std::chrono::steady_clock::now();

    std::cout << "boost::charconv::to_chars<" << boost::core::type_name<T>() << ">, " << label << ", " << precision << ": " << std::setw( 5 ) << ( t2 - t1 ) / 1ms << " ms (s=" << s << ")\n";
}

template<class T> static void test()
{
    std::vector<T> data;
    init_input_data( data );

    test_snprintf( data, false, "scientific", 0 );
    test_std_to_chars( data, false, "scientific", 0 );
    test_boost_to_chars( data, false, "scientific", 0 );

    std::cout << std::endl;

    test_snprintf( data, false, "scientific", 6 );
    test_std_to_chars( data, false, "scientific", 6 );
    test_boost_to_chars( data, false, "scientific", 6 );

    std::cout << std::endl;

    test_snprintf( data, true, "general", 0 );
    test_std_to_chars( data, true, "general", 0 );
    test_boost_to_chars( data, true, "general", 0 );

    std::cout << std::endl;

    test_snprintf( data, true, "general", 6 );
    test_std_to_chars( data, true, "general", 6 );
    test_boost_to_chars( data, true, "general", 6 );

    std::cout << std::endl;
}

#ifdef BOOST_CHARCONV_HAS_STDFLOAT128
template<> void test<std::float128_t>()
{
    std::vector<std::float128_t> data;
    init_input_data( data );

    test_std_to_chars( data, false, "scientific", 0 );
    test_boost_to_chars( data, false, "scientific", 0 );

    std::cout << std::endl;

    test_std_to_chars( data, false, "scientific", 6 );
    test_boost_to_chars( data, false, "scientific", 6 );

    std::cout << std::endl;

    test_std_to_chars( data, true, "general", 0 );
    test_boost_to_chars( data, true, "general", 0 );

    std::cout << std::endl;

    test_std_to_chars( data, true, "general", 6 );
    test_boost_to_chars( data, true, "general", 6 );

    std::cout << std::endl;
}
#endif

int main()
{
    std::cout << BOOST_COMPILER << "\n";
    std::cout << BOOST_STDLIB << "\n\n";

    test<float>();
    test<double>();
    #ifdef BOOST_CHARCONV_HAS_STDFLOAT128
    test<std::float128_t>();
    #endif
}
