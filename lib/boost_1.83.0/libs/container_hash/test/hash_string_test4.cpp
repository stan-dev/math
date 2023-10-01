// Copyright 2021, 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/container_hash/hash.hpp>
#include <boost/core/lightweight_test.hpp>
#include <string>
#include <cstddef>

void test( std::size_t n, unsigned char ch )
{
    typedef boost::hash<std::string> hash;

    std::string const v( n, static_cast<char>( ch ) );

    for( std::size_t i = 0; i < n * 8; ++i )
    {
        std::string w( v );

        unsigned char ch2 = static_cast<unsigned char>( w[ i / 8 ] );

        ch2 = static_cast<unsigned char>( ch2 ^ ( 1 << ( i % 8 ) ) );

        w[ i / 8 ] = static_cast<char>( ch2 );

        BOOST_TEST_NE( hash()( v ), hash()( w ) );
    }
}

int main()
{
    for( unsigned ch = 0; ch < 256; ++ch )
    {
        for( std::size_t n = 1; n < 32; ++n )
        {
            test( n, static_cast<unsigned char>( ch ) );
        }
    }

    return boost::report_errors();
}
