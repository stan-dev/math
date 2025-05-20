// Copyright 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt

#include <boost/system/error_code.hpp>
#include <boost/system/generic_category.hpp>
#include <boost/system/system_category.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cerrno>
#include <system_error>

namespace sys = boost::system;

int main()
{
    // normal against normal (equal, system)
    {
        sys::error_code e2( 0, sys::system_category() );
        sys::error_code e3( e2.value(), e2.category() );

        BOOST_TEST( e2 != e3 || hash_value( e2 ) == hash_value( e3 ) );
        BOOST_TEST( e2 == e3 || hash_value( e2 ) != hash_value( e3 ) );
    }

    // normal against normal (equal, generic)
    {
        sys::error_code e2( EINVAL, sys::generic_category() );
        sys::error_code e3( e2.value(), e2.category() );

        BOOST_TEST( e2 != e3 || hash_value( e2 ) == hash_value( e3 ) );
        BOOST_TEST( e2 == e3 || hash_value( e2 ) != hash_value( e3 ) );
    }

    // normal against normal (inequal, value, generic)
    {
        sys::error_code e2( 0, sys::generic_category() );
        sys::error_code e3( EINVAL, sys::generic_category() );

        BOOST_TEST( e2 != e3 || hash_value( e2 ) == hash_value( e3 ) );
        BOOST_TEST( e2 == e3 || hash_value( e2 ) != hash_value( e3 ) );
    }

    // normal against normal (inequal, value, system)
    {
        sys::error_code e2( 1, sys::system_category() );
        sys::error_code e3( 2, sys::system_category() );

        BOOST_TEST( e2 != e3 || hash_value( e2 ) == hash_value( e3 ) );
        BOOST_TEST( e2 == e3 || hash_value( e2 ) != hash_value( e3 ) );
    }

    // normal against normal (inequal, category)
    {
        sys::error_code e2( 0, sys::system_category() );
        sys::error_code e3( 0, sys::generic_category() );

        BOOST_TEST( e2 != e3 || hash_value( e2 ) == hash_value( e3 ) );
        BOOST_TEST( e2 == e3 || hash_value( e2 ) != hash_value( e3 ) );
    }

    // empty against normal
    {
        sys::error_code e2;
        sys::error_code e3( e2.value(), e2.category() );

        BOOST_TEST( e2 != e3 || hash_value( e2 ) == hash_value( e3 ) );
        BOOST_TEST( e2 == e3 || hash_value( e2 ) != hash_value( e3 ) );
    }

    // std:: wrapping against normal
    {
        std::error_code e1( EINVAL, std::generic_category() );

        sys::error_code e2( e1 );
        sys::error_code e3( e2.value(), e2.category() );

        BOOST_TEST( e2 != e3 || hash_value( e2 ) == hash_value( e3 ) );
        BOOST_TEST( e2 == e3 || hash_value( e2 ) != hash_value( e3 ) );
    }

    // empty against wrapping std:: empty
    {
        std::error_code e1;

        sys::error_code e2( e1 );
        sys::error_code e3;

        BOOST_TEST( e2 != e3 || hash_value( e2 ) == hash_value( e3 ) );
        BOOST_TEST( e2 == e3 || hash_value( e2 ) != hash_value( e3 ) );
    }

    // empty against roundtrip via std
    {
        sys::error_code e2;

        std::error_code e1( e2 );
        sys::error_code e3( e1 );

        BOOST_TEST( e2 != e3 || hash_value( e2 ) == hash_value( e3 ) );
        BOOST_TEST( e2 == e3 || hash_value( e2 ) != hash_value( e3 ) );
    }

    // normal/generic against roundtrip via std
    {
        sys::error_code e2( EINVAL, boost::system::generic_category() );

        std::error_code e1( e2 );
        sys::error_code e3( e1 );

        BOOST_TEST( e2 != e3 || hash_value( e2 ) == hash_value( e3 ) );
        BOOST_TEST( e2 == e3 || hash_value( e2 ) != hash_value( e3 ) );
    }

    // normal/system against roundtrip via std
    {
        sys::error_code e2( 0, boost::system::system_category() );

        std::error_code e1( e2 );
        sys::error_code e3( e1 );

        BOOST_TEST( e2 != e3 || hash_value( e2 ) == hash_value( e3 ) );
        BOOST_TEST( e2 == e3 || hash_value( e2 ) != hash_value( e3 ) );
    }

    return boost::report_errors();
}
