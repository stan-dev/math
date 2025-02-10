// Copyright 2023 Peter Dimov
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

constexpr unsigned N = 2'000'000;
constexpr int K = 25;

template<class T> static BOOST_NOINLINE void init_input_data( std::vector<std::string>& data )
{
    data.reserve( N );

    boost::detail::splitmix64 rng;

    for( unsigned i = 0; i < N; ++i )
    {
        T x = static_cast<T>( static_cast<typename std::make_unsigned<T>::type>( rng() ) );

        char buffer[ 21 ];
        auto r = boost::charconv::to_chars( buffer, buffer + sizeof( buffer ), x );

        std::string y( buffer, r.ptr );
        data.push_back( y );
    }
}

using namespace std::chrono_literals;

template<class T> static void BOOST_NOINLINE test_std_from_chars( std::vector<std::string> const& data )
{
    auto t1 = std::chrono::steady_clock::now();

    std::size_t s = 0;

    for( int i = 0; i < K; ++i )
    {
        for( auto const& x: data )
        {
            T y;
            std::from_chars( x.data(), x.data() + x.size(), y );

            s += static_cast<std::size_t>( y );
        }
    }

    auto t2 = std::chrono::steady_clock::now();

    std::cout << "            std::from_chars<" << boost::core::type_name<T>() << ">: " << std::setw( 5 ) << ( t2 - t1 ) / 1ms << " ms (s=" << s << ")\n";
}

template<class T> static void BOOST_NOINLINE test_boost_from_chars( std::vector<std::string> const& data )
{
    auto t1 = std::chrono::steady_clock::now();

    std::size_t s = 0;

    for( int i = 0; i < K; ++i )
    {
        for( auto const& x: data )
        {
            T y;
            boost::charconv::from_chars( x.data(), x.data() + x.size(), y );

            s += static_cast<std::size_t>( y );
        }
    }

    auto t2 = std::chrono::steady_clock::now();

    std::cout << "boost::charconv::from_chars<" << boost::core::type_name<T>() << ">: " << std::setw( 5 ) << ( t2 - t1 ) / 1ms << " ms (s=" << s << ")\n";
}

template<class T> static void test()
{
    std::vector<std::string> data;
    init_input_data<T>( data );

    test_std_from_chars<T>( data );
    test_boost_from_chars<T>( data );

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
