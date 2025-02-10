// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

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
        data.push_back( static_cast<T>( static_cast<typename std::make_unsigned<T>::type>( rng() ) ) );
    }
}

using namespace std::chrono_literals;

template<class T> static BOOST_NOINLINE void test_snprintf( std::vector<T> const& data )
{
    auto t1 = std::chrono::steady_clock::now();

    std::size_t s = 0;

    char const* format = "%d";

    BOOST_IF_CONSTEXPR( std::is_same<T, unsigned>::value )
    {
        format = "%u";
    }

    BOOST_IF_CONSTEXPR( std::is_same<T, long>::value )
    {
        format = "%ld";
    }

    BOOST_IF_CONSTEXPR( std::is_same<T, unsigned long>::value )
    {
        format = "%lu";
    }

    BOOST_IF_CONSTEXPR( std::is_same<T, long long>::value )
    {
        format = "%lld";
    }

    BOOST_IF_CONSTEXPR( std::is_same<T, unsigned long long>::value )
    {
        format = "%llu";
    }

    for( int i = 0; i < K; ++i )
    {
        char buffer[ 22 ];

        for( auto x: data )
        {
            auto r = std::snprintf( buffer, sizeof( buffer ), format, x );
            s += r;
            s += static_cast<unsigned char>( buffer[0] );
        }
    }

    auto t2 = std::chrono::steady_clock::now();

    std::cout << "            std::snprintf<" << boost::core::type_name<T>() << ">: " << std::setw( 5 ) << ( t2 - t1 ) / 1ms << " ms (s=" << s << ")\n";
}

template<class T> static BOOST_NOINLINE void test_std_to_chars( std::vector<T> const& data )
{
    auto t1 = std::chrono::steady_clock::now();

    std::size_t s = 0;

    for( int i = 0; i < K; ++i )
    {
        char buffer[ 21 ];

        for( auto x: data )
        {
            auto r = std::to_chars( buffer, buffer + sizeof( buffer ), x );
            s += static_cast<std::size_t>( r.ptr - buffer );
            s += static_cast<unsigned char>( buffer[0] );
        }
    }

    auto t2 = std::chrono::steady_clock::now();

    std::cout << "            std::to_chars<" << boost::core::type_name<T>() << ">: " << std::setw( 5 ) << ( t2 - t1 ) / 1ms << " ms (s=" << s << ")\n";
}

template<class T> static BOOST_NOINLINE void test_boost_to_chars( std::vector<T> const& data )
{
    auto t1 = std::chrono::steady_clock::now();

    std::size_t s = 0;

    for( int i = 0; i < K; ++i )
    {
        char buffer[ 21 ];

        for( auto x: data )
        {
            auto r = boost::charconv::to_chars( buffer, buffer + sizeof( buffer ), x );
            s += static_cast<std::size_t>( r.ptr - buffer );
            s += static_cast<unsigned char>( buffer[0] );
        }
    }

    auto t2 = std::chrono::steady_clock::now();

    std::cout << "boost::charconv::to_chars<" << boost::core::type_name<T>() << ">: " << std::setw( 5 ) << ( t2 - t1 ) / 1ms << " ms (s=" << s << ")\n";
}

template<class T> static void test()
{
    std::vector<T> data;
    init_input_data( data );

    test_snprintf( data );
    test_std_to_chars( data );
    test_boost_to_chars( data );

    std::cout << std::endl;
}

int main()
{
    std::cout << BOOST_COMPILER << "\n";
    std::cout << BOOST_STDLIB << "\n\n";

    test<signed char>();
    test<unsigned char>();

    test<short>();
    test<unsigned short>();

    test<int>();
    test<unsigned int>();

    test<long>();
    test<unsigned long>();

    test<long long>();
    test<unsigned long long>();
}
