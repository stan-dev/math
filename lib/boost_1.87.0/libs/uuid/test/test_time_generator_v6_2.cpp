// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/time_generator_v6.hpp>
#include <boost/uuid/uuid_clock.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/config.hpp>
#include <chrono>
#include <atomic>
#include <set>

#if defined(BOOST_LIBSTDCXX_VERSION) && BOOST_LIBSTDCXX_VERSION >= 130000
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=114865
# pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif

using namespace boost::uuids;

uuid generate_and_test( time_generator_v6& gen )
{
    auto sys_before = std::chrono::time_point_cast<uuid_clock::duration>( std::chrono::system_clock::now() );

    uuid u = gen();

    BOOST_TEST_EQ( u.variant(), uuid::variant_rfc_4122 );
    BOOST_TEST_EQ( u.version(), uuid::version_time_based_v6 );

    auto sys_after = std::chrono::time_point_cast<uuid_clock::duration>( std::chrono::system_clock::now() );

    auto uuid_time_point = u.time_point_v6();

    auto sys_time_point = uuid_clock::to_sys( uuid_time_point );

    BOOST_TEST( sys_before <= sys_time_point );
    BOOST_TEST( sys_time_point <= sys_after );

    return u;
}

int main()
{
    int const N = 1024;

    {
        std::set<uuid> set;

        uuid::node_type node{{ 0x01, 0x02, 0x03, 0x04, 0x05, 0x06 }};
        std::atomic< time_generator_v6::state_type > state{{ 0, 0x2222 }};

        time_generator_v6 gen( node, state );

        uuid u1 = generate_and_test( gen );

        BOOST_TEST( u1.node_identifier() == node );
        BOOST_TEST_EQ( u1.clock_seq(), 0x2222 );

        set.insert( u1 );

        for( int i = 0; i < N; ++i )
        {
            uuid u2 = generate_and_test( gen );

            BOOST_TEST( u2.node_identifier() == node );

            set.insert( u2 );
        }

        BOOST_TEST_EQ( set.size(), N + 1 );
    }

    return boost::report_errors();
}
