// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/time_generator_v6.hpp>
#include <boost/uuid/uuid_clock.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/config.hpp>
#include <atomic>
#include <chrono>
#include <thread>
#include <set>
#include <vector>

#if defined(BOOST_LIBSTDCXX_VERSION) && BOOST_LIBSTDCXX_VERSION >= 130000
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=114865
# pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif

using namespace boost::uuids;

uuid::node_type const node{{ 0x01, 0x02, 0x03, 0x04, 0x05, 0x06 }};

int const N = 1024;

void threadfunc( std::atomic<time_generator_v6::state_type>& state, std::vector<uuid>& v )
{
    time_generator_v6 gen( node, state );

    for( int i = 0; i < N; ++i )
    {
        v.push_back( gen() );
    }
}

int main()
{
    int const M = 8;

    std::thread th[ M ];
    std::vector<uuid> v[ M ];

    std::atomic<time_generator_v6::state_type> state{{ 0, 0 }};

    for( int i = 0; i < M; ++i )
    {
        th[ i ] = std::thread( threadfunc, std::ref( state ), std::ref( v[ i ] ) );
    }

    for( int i = 0; i < M; ++i )
    {
        th[ i ].join();
    }

    std::set<uuid> set;

    for( int i = 0; i < M; ++i )
    {
        set.insert( v[ i ].begin(), v[ i ].end() );
    }

    BOOST_TEST_EQ( set.size(), N * M );

    return boost::report_errors();
}
