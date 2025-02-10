// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/string_generator.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::uuids;

int main()
{
    // Test vectors from IETF RFC 9562, Appendix A
    // https://datatracker.ietf.org/doc/rfc9562/

    string_generator uuid_from_string;

    // A.1, UUIDv1
    {
        uuid::node_type const node{{ 0x9F, 0x6B, 0xDE, 0xCE, 0xD8, 0x46 }};

        uuid u = uuid_from_string( "C232AB00-9414-11EC-B3C8-9F6BDECED846" );

        BOOST_TEST_EQ( u.variant(), uuid::variant_rfc_4122 );
        BOOST_TEST_EQ( u.version(), uuid::version_time_based );
        BOOST_TEST_EQ( u.timestamp_v1(), 0x1EC9414C232AB00 );
        BOOST_TEST_EQ( u.clock_seq(), 0x33C8 );
        BOOST_TEST( u.node_identifier() == node );
    }

    // A.2, UUIDv3
    {
        uuid u = uuid_from_string( "5df41881-3aed-3515-88a7-2f4a814cf09e" );

        BOOST_TEST_EQ( u.variant(), uuid::variant_rfc_4122 );
        BOOST_TEST_EQ( u.version(), uuid::version_name_based_md5 );
    }

    // A.3, UUIDv4
    {
        uuid u = uuid_from_string( "919108f7-52d1-4320-9bac-f847db4148a8" );

        BOOST_TEST_EQ( u.variant(), uuid::variant_rfc_4122 );
        BOOST_TEST_EQ( u.version(), uuid::version_random_number_based );
    }

    // A.4, UUIDv5
    {
        uuid u = uuid_from_string( "2ed6657d-e927-568b-95e1-2665a8aea6a2" );

        BOOST_TEST_EQ( u.variant(), uuid::variant_rfc_4122 );
        BOOST_TEST_EQ( u.version(), uuid::version_name_based_sha1 );
    }

    // A.5, UUIDv6
    {
        uuid::node_type const node{{ 0x9F, 0x6B, 0xDE, 0xCE, 0xD8, 0x46 }};

        uuid u = uuid_from_string( "1EC9414C-232A-6B00-B3C8-9F6BDECED846" );

        BOOST_TEST_EQ( u.variant(), uuid::variant_rfc_4122 );
        BOOST_TEST_EQ( u.version(), uuid::version_time_based_v6 );
        BOOST_TEST_EQ( u.timestamp_v6(), 0x1EC9414C232AB00 );
        BOOST_TEST_EQ( u.clock_seq(), 0x33C8 );
        BOOST_TEST( u.node_identifier() == node );
    }

    // A.6, UUIDv7
    {
        uuid u = uuid_from_string( "017F22E2-79B0-7CC3-98C4-DC0C0C07398F" );

        BOOST_TEST_EQ( u.variant(), uuid::variant_rfc_4122 );
        BOOST_TEST_EQ( u.version(), uuid::version_time_based_v7 );
        BOOST_TEST_EQ( u.timestamp_v7(), 0x017F22E279B0 );
    }

    // B.1, UUIDv8
    {
        uuid u = uuid_from_string( "2489E9AD-2EE2-8E00-8EC9-32D5F69181C0" );

        BOOST_TEST_EQ( u.variant(), uuid::variant_rfc_4122 );
        BOOST_TEST_EQ( u.version(), uuid::version_custom_v8 );
    }

    // B.2, UUIDv8
    {
        uuid u = uuid_from_string( "5c146b14-3c52-8afd-938a-375d0df1fbf6" );

        BOOST_TEST_EQ( u.variant(), uuid::variant_rfc_4122 );
        BOOST_TEST_EQ( u.version(), uuid::version_custom_v8 );
    }

    return boost::report_errors();
}
