// Copyright 2021, 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/container_hash/hash.hpp>
#include <boost/core/lightweight_test.hpp>
#include <string>

void test( unsigned char ch )
{
    typedef boost::hash<std::string> hash;

    int const N = 32;

    std::size_t h[ N ];

    std::string v;

    for( int i = 0; i < N; ++i )
    {
        h[ i ] = hash()( v );

        BOOST_TEST_EQ( h[ i ], hash()( v ) );

        for( int j = 0; j < i; ++j )
        {
            BOOST_TEST_NE( h[ j ], h[ i ] );
        }

        v.push_back( static_cast<char>( ch ) );
    }
}

int main()
{
    for( unsigned ch = 0; ch < 256; ++ch )
    {
        test( static_cast<unsigned char>( ch ) );
    }

    return boost::report_errors();
}
