// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/uuid.hpp>
#include <boost/core/lightweight_test.hpp>
#include <unordered_set>
#include <algorithm>
#include <cstdio>

int main()
{
    using namespace boost::uuids;

    for( unsigned ch = 0; ch < 256; ++ch )
    {
        std::unordered_set<std::size_t> set;

        uuid u1 = {};

        std::fill( u1.data + 0, u1.data + 16, static_cast<unsigned char>( ch ) );

        set.insert( hash_value( u1 ) );

        for( int i = 0; i < 16; ++i )
        {
            for( int j = 0; j < 8; ++j )
            {
                uuid u2 = u1;

                u2.data[ i ] ^= static_cast<unsigned char>( 1 << j );

                set.insert( hash_value( u2 ) );
            }
        }

        if( !BOOST_TEST_EQ( set.size(), 16 * 8 + 1 ) )
        {
            std::printf( "^-- Collisions detected with ch=0x%02X\n\n", ch );
        }
    }

    return boost::report_errors();
}
