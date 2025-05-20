// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/time_generator_v7.hpp>
#include <boost/core/lightweight_test.hpp>
#include <chrono>
#include <set>

using namespace boost::uuids;

uuid generate_and_test( time_generator_v7& gen )
{
    auto sys_before = std::chrono::time_point_cast<std::chrono::milliseconds>( std::chrono::system_clock::now() );

    uuid u = gen();

    BOOST_TEST_EQ( u.variant(), uuid::variant_rfc_4122 );
    BOOST_TEST_EQ( u.version(), uuid::version_time_based_v7 );

    // The below is an unsupported way to obtain the
    // microsecond part of our UUIDv7 timestamps.
    //
    // It's not guaranteed by the specification and may
    // break in future releases.

    std::uint64_t time_in_us = u.timestamp_v6() & 0xFFF; // v6 time_low field

    // time_in_us occasionally exceeds 999 when the clock
    // resolution is sufficiently low (MinGW)

    BOOST_TEST_LT( time_in_us, 1024 );

    auto sys_after = std::chrono::time_point_cast<std::chrono::milliseconds>( std::chrono::system_clock::now() );

    auto sys_time_point = u.time_point_v7();

    BOOST_TEST( sys_before <= sys_time_point );
    BOOST_TEST( sys_time_point <= sys_after );

    return u;
}

int main()
{
    int const N = 1024;

    {
        std::set<uuid> set;

        time_generator_v7 gen;

        uuid u1 = generate_and_test( gen );

        set.insert( u1 );

        for( int i = 0; i < N; ++i )
        {
            uuid u2 = generate_and_test( gen );

            set.insert( u2 );
        }

        BOOST_TEST_EQ( set.size(), N + 1 );
    }

    return boost::report_errors();
}
