// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/time_generator_v7.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/core/lightweight_test.hpp>
#include <utility>

// Test that the UUIDs produced by copies of time_generator_v7
// aren't duplicates

using namespace boost::uuids;

int main()
{
    int const N = 16;

    for( int i = 0; i < N; ++i )
    {
        time_generator_v7 gen;
        time_generator_v7 gen2( gen );

        uuid u1 = gen();
        uuid u2 = gen2();

        BOOST_TEST_NE( u1, u2 );
    }

    for( int i = 0; i < N; ++i )
    {
        time_generator_v7 gen;
        time_generator_v7 gen2( std::move( gen ) );

        uuid u1 = gen();
        uuid u2 = gen2();

        BOOST_TEST_NE( u1, u2 );
    }

    for( int i = 0; i < N; ++i )
    {
        time_generator_v7 gen;
        time_generator_v7 gen2;

        gen2 = gen;

        uuid u1 = gen();
        uuid u2 = gen2();

        BOOST_TEST_NE( u1, u2 );
    }

    for( int i = 0; i < N; ++i )
    {
        time_generator_v7 gen;
        time_generator_v7 gen2;

        gen2 = std::move( gen );

        uuid u1 = gen();
        uuid u2 = gen2();

        BOOST_TEST_NE( u1, u2 );
    }

    return boost::report_errors();
}
