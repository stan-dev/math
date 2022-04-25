// Copyright 2021, 2022 Peter Dimov
// Distributed under the Boost Software License, Version 1.0
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/system/system_error.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cerrno>

namespace sys = boost::system;

int main()
{
    {
        sys::error_code ec( 5, sys::generic_category() );
        sys::system_error x1( ec );

        BOOST_TEST_EQ( std::string( x1.what() ), ec.what() );
    }

    {
        sys::error_code ec( 5, sys::generic_category() );
        sys::system_error x1( ec, "prefix" );

        BOOST_TEST_EQ( std::string( x1.what() ), "prefix: " + ec.what() );
    }

    {
        sys::error_code ec( 5, sys::generic_category() );
        sys::system_error x1( ec, std::string( "prefix2" ) );

        BOOST_TEST_EQ( std::string( x1.what() ), "prefix2: " + ec.what() );
    }

    {
        BOOST_STATIC_CONSTEXPR boost::source_location loc = BOOST_CURRENT_LOCATION;

        sys::error_code ec( 5, sys::generic_category(), &loc );
        sys::system_error x1( ec, "prefix3" );

        BOOST_TEST_EQ( std::string( x1.what() ), "prefix3: " + ec.what() );
    }

    {
        sys::system_error x1( 5, sys::system_category() );

        BOOST_TEST_EQ( std::string( x1.what() ), sys::error_code( 5, sys::system_category() ).what() );
    }

    {
        sys::system_error x1( 5, sys::system_category(), "prefix" );

        BOOST_TEST_EQ( std::string( x1.what() ), "prefix: " + sys::error_code( 5, sys::system_category() ).what() );
    }

    {
        sys::system_error x1( 5, sys::system_category(), std::string( "prefix2" ) );

        BOOST_TEST_EQ( std::string( x1.what() ), "prefix2: " + sys::error_code( 5, sys::system_category() ).what() );
    }

    return boost::report_errors();
}
