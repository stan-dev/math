// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/time_generator.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::uuids;

int main()
{
    uuid u1 = time_generator_v1()();

    BOOST_TEST_EQ( u1.variant(), uuid::variant_rfc_4122 );
    BOOST_TEST_EQ( u1.version(), uuid::version_time_based );

    uuid u2 = time_generator_v6()();

    BOOST_TEST_EQ( u2.variant(), uuid::variant_rfc_4122 );
    BOOST_TEST_EQ( u2.version(), uuid::version_time_based_v6 );

    uuid u3 = time_generator_v7()();

    BOOST_TEST_EQ( u3.variant(), uuid::variant_rfc_4122 );
    BOOST_TEST_EQ( u3.version(), uuid::version_time_based_v7 );

    return boost::report_errors();
}
