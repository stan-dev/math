// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/time_generator_v7.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/core/lightweight_test.hpp>
#include <vector>

// Test that the UUIDs produced by time_generator_v7
// are strictly monotonic

using namespace boost::uuids;

int main()
{
    int const N = 1024;

    {
        std::vector<uuid> v( N );

        time_generator_v7 gen;

        for( int i = 0; i < N; ++i )
        {
            v[ i ] = gen();
        }

        for( int i = 0; i < N - 1; ++i )
        {
            BOOST_TEST_LT( v[ i ], v[ i+1 ] );
        }
    }

    return boost::report_errors();
}
