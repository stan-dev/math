// Copyright 2023 Peter Dimov
// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/charconv/from_chars.hpp>
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
#include <random>

constexpr unsigned N = 2'000'000;
constexpr int K = 10;

template<class T> static BOOST_NOINLINE void init_input_data( std::vector<std::string>& data, bool general )
{
    data.reserve( N );

    boost::detail::splitmix64 rng;

    for( unsigned i = 0; i < N; ++i )
    {
        std::uint64_t tmp = rng();

        T x;
        std::memcpy( &x, &tmp, sizeof(x) );

        if( !std::isfinite(x) ) continue;

        char buffer[ 64 ];
        auto r = boost::charconv::to_chars( buffer, buffer + sizeof( buffer ), x, general? boost::charconv::chars_format::general: boost::charconv::chars_format::scientific );

        std::string y( buffer, r.ptr );
        data.push_back( y );
    }
}

template <>
BOOST_NOINLINE void init_input_data<long double>( std::vector<std::string>& data, bool general )
{
    data.reserve(N / 10);

    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_real_distribution<long double> dist(0.0L, (std::numeric_limits<long double>::max)());

    for( unsigned i = 0; i < N / 10; ++i )
    {
        const long double x = dist(rng);

        char buffer[ 64 ];
        auto r = boost::charconv::to_chars( buffer, buffer + sizeof( buffer ), x, general? boost::charconv::chars_format::general : boost::charconv::chars_format::scientific );

        std::string y( buffer, r.ptr );
        data.push_back( y );
    }
}

template<class T> static BOOST_NOINLINE void init_input_data_uint64( std::vector<std::string>& data )
{
    data.reserve( N );

    boost::detail::splitmix64 rng;

    for( unsigned i = 0; i < N; ++i )
    {
        std::uint64_t x = rng();

        char buffer[ 64 ];
        auto r = boost::charconv::to_chars( buffer, buffer + sizeof( buffer ), x );

        std::string y( buffer, r.ptr );
        data.push_back( y );
    }
}

using namespace std::chrono_literals;

template<class T> void test_strtox( std::vector<std::string> const& data, bool general, char const* label );

template<> BOOST_NOINLINE void test_strtox<float>( std::vector<std::string> const& data, bool, char const* label )
{
    auto t1 = std::chrono::steady_clock::now();

    double s = 0;

    for( int i = 0; i < K; ++i )
    {
        for( auto const& x: data )
        {
            float y = std::strtof( x.c_str(), nullptr );
            s = s / 16.0 + y;
        }
    }

    auto t2 = std::chrono::steady_clock::now();

    std::cout << "                std::strtox<float>, " << label << ": " << std::setw( 5 ) << ( t2 - t1 ) / 1ms << " ms (s=" << s << ")\n";
}

template<> BOOST_NOINLINE void test_strtox<double>( std::vector<std::string> const& data, bool, char const* label )
{
    auto t1 = std::chrono::steady_clock::now();

    double s = 0;

    for( int i = 0; i < K; ++i )
    {
        for( auto const& x: data )
        {
            double y = std::strtod( x.c_str(), nullptr );
            s = s / 16.0 + y;
        }
    }

    auto t2 = std::chrono::steady_clock::now();

    std::cout << "                std::strtox<double>, " << label << ": " << std::setw( 5 ) << ( t2 - t1 ) / 1ms << " ms (s=" << s << ")\n";
}

template<> BOOST_NOINLINE void test_strtox<long double>( std::vector<std::string> const& data, bool, char const* label )
{
    auto t1 = std::chrono::steady_clock::now();

    long double s = 0;

    for( int i = 0; i < K; ++i )
    {
        for( auto const& x: data )
        {
            double y = std::strtold( x.c_str(), nullptr );
            s = s / 16.0L + y;
        }
    }

    auto t2 = std::chrono::steady_clock::now();

    std::cout << "                std::strtox<long double>, " << label << ": " << std::setw( 5 ) << ( t2 - t1 ) / 1ms << " ms (s=" << s << ")\n";
}

template<class T> static BOOST_NOINLINE void test_std_from_chars( std::vector<std::string> const& data, bool general, char const* label )
{
    auto t1 = std::chrono::steady_clock::now();

    double s = 0;

    for( int i = 0; i < K; ++i )
    {
        for( auto const& x: data )
        {
            T y;
            std::from_chars( x.data(), x.data() + x.size(), y, general? std::chars_format::general: std::chars_format::scientific );

            s = s / 16.0 + y;
        }
    }

    auto t2 = std::chrono::steady_clock::now();

    std::cout << "            std::from_chars<" << boost::core::type_name<T>() << ">, " << label << ": " << std::setw( 5 ) << ( t2 - t1 ) / 1ms << " ms (s=" << s << ")\n";
}

template<> BOOST_NOINLINE void test_std_from_chars<long double>( std::vector<std::string> const& data, bool general, char const* label )
{
    auto t1 = std::chrono::steady_clock::now();

    long double s = 0;

    for( int i = 0; i < K; ++i )
    {
        for( auto const& x: data )
        {
            long double y;
            std::from_chars( x.data(), x.data() + x.size(), y, general? std::chars_format::general: std::chars_format::scientific );

            s = s / 16.0L + y;
        }
    }

    auto t2 = std::chrono::steady_clock::now();

    std::cout << "            std::from_chars<long double>, " << label << ": " << std::setw( 5 ) << ( t2 - t1 ) / 1ms << " ms (s=" << s << ")\n";
}

template<class T> static BOOST_NOINLINE void test_boost_from_chars( std::vector<std::string> const& data, bool general, char const* label )
{
    auto t1 = std::chrono::steady_clock::now();

    double s = 0;

    for( int i = 0; i < K; ++i )
    {
        for( auto const& x: data )
        {
            T y;
            boost::charconv::from_chars( x.data(), x.data() + x.size(), y, general? boost::charconv::chars_format::general: boost::charconv::chars_format::scientific );

            s = s / 16.0 + y;
        }
    }

    auto t2 = std::chrono::steady_clock::now();

    std::cout << "boost::charconv::from_chars<" << boost::core::type_name<T>() << ">, " << label << ": " << std::setw( 5 ) << ( t2 - t1 ) / 1ms << " ms (s=" << s << ")\n";
}

template<>
BOOST_NOINLINE void test_boost_from_chars<long double>( std::vector<std::string> const& data, bool general, char const* label )
{
    auto t1 = std::chrono::steady_clock::now();

    long double s = 0;

    for( int i = 0; i < K; ++i )
    {
        for( auto const& x: data )
        {
            long double y;
            boost::charconv::from_chars( x.data(), x.data() + x.size(), y, general? boost::charconv::chars_format::general: boost::charconv::chars_format::scientific );

            s = s / 16.0L + y;
        }
    }

    auto t2 = std::chrono::steady_clock::now();

    std::cout << "boost::charconv::from_chars<long double>, " << label << ": " << std::setw( 5 ) << ( t2 - t1 ) / 1ms << " ms (s=" << s << ")\n";
}

template<class T> static void test( bool general )
{
    std::vector<std::string> data;
    init_input_data<T>( data, general );

    char const* label = general? "general": "scientific";

    test_strtox<T>( data, general, label );
    test_std_from_chars<T>( data, general, label );
    test_boost_from_chars<T>( data, general, label );

    std::cout << std::endl;
}

template<class T> static void test2()
{
    std::vector<std::string> data;
    init_input_data_uint64<T>( data );

    bool general = true;
    char const* label = "uint64";

    test_strtox<T>( data, general, label );
    test_std_from_chars<T>( data, general, label );
    test_boost_from_chars<T>( data, general, label );

    std::cout << std::endl;
}

int main()
{
    std::cout << BOOST_COMPILER << "\n";
    std::cout << BOOST_STDLIB << "\n\n";

    test<float>( false );
    test<double>( false );
    test<long double>( false );

    test<float>( true );
    test<double>( true );
    test<long double>( true );

    test2<float>();
    test2<double>();
}
