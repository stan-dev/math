// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/unordered/detail/foa/rw_spinlock.hpp>
#include <boost/core/lightweight_test.hpp>
#include <mutex>
#include <thread>
#include <cstdio>

using boost::unordered::detail::foa::rw_spinlock;

static int count = 0;
static rw_spinlock sp;

void f( int k, int n )
{
    std::printf( "Thread %d started.\n", k );

    for( int i = 0; i < n; ++i )
    {
        std::lock_guard<rw_spinlock> lock( sp );
        ++count;
    }

    std::printf( "Thread %d finished.\n", k );
}

int main()
{
    int const N = 1000000; // iterations
    int const M = 8;       // threads

    std::thread th[ M ];

    for( int i = 0; i < M; ++i )
    {
        th[ i ] = std::thread( f, i, N );
    }

    for( int i = 0; i < M; ++i )
    {
        th[ i ].join();
    }

    BOOST_TEST_EQ( count, N * M );

    return boost::report_errors();
}
