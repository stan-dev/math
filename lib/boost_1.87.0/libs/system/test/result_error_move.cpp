// Copyright 2023 Klemens Morgenstern
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/system/result.hpp>
#include <boost/core/lightweight_test.hpp>
#include <string>

using namespace boost::system;

struct X
{
    int v_;

    explicit X( int v = 0 ): v_( v ) {}

    X( X const& ) = delete;
    X& operator=( X const& ) = delete;

    X( X && ) = default;
    X& operator=( X && ) = default;
};

int main()
{
    {
        result<std::string, X> r( 1 );

        BOOST_TEST( !r.has_value() );
        BOOST_TEST( r.has_error() );

        BOOST_TEST_EQ( std::move( r ).error().v_, 1 );
    }

    {
        BOOST_TEST(( !result<std::string, X>( 1 ).has_value() ));
        BOOST_TEST(( result<std::string, X>( 1 ).has_error() ));

        BOOST_TEST_EQ( (result<std::string, X>( 1 ).error().v_), 1 );
    }

    {
        result<std::string, X> r( "s" );

        BOOST_TEST( r.has_value() );
        BOOST_TEST( !r.has_error() );

        BOOST_TEST_EQ( std::move( r ).error().v_, 0 );
    }

    {
        BOOST_TEST(( result<std::string, X>( "s" ).has_value() ));
        BOOST_TEST(( !result<std::string, X>( "s" ).has_error() ));

        BOOST_TEST_EQ( (result<std::string, X>( "s" ).error().v_), 0 );
    }

    //

    {
        result<void, X> r( 1 );

        BOOST_TEST( !r.has_value() );
        BOOST_TEST( r.has_error() );

        BOOST_TEST_EQ( std::move( r ).error().v_, 1 );
    }

    {
        BOOST_TEST(( !result<void, X>( 1 ).has_value() ));
        BOOST_TEST(( result<void, X>( 1 ).has_error() ));

        BOOST_TEST_EQ( (result<void, X>( 1 ).error().v_), 1 );
    }

    {
        result<void, X> r;

        BOOST_TEST( r.has_value() );
        BOOST_TEST( !r.has_error() );

        BOOST_TEST_EQ( std::move( r ).error().v_, 0 );
    }

    {
        BOOST_TEST(( result<void, X>().has_value() ));
        BOOST_TEST(( !result<void, X>().has_error() ));

        BOOST_TEST_EQ( (result<void, X>().error().v_), 0 );
    }

    //

    {
        result<double&, X> r( 1 );

        BOOST_TEST( !r.has_value() );
        BOOST_TEST( r.has_error() );

        BOOST_TEST_EQ( std::move( r ).error().v_, 1 );
    }

    {
        BOOST_TEST(( !result<double&, X>( 1 ).has_value() ));
        BOOST_TEST(( result<double&, X>( 1 ).has_error() ));

        BOOST_TEST_EQ( (result<double&, X>( 1 ).error().v_), 1 );
    }

    {
        double x = 1.0;

        result<double&, X> r( x );

        BOOST_TEST( r.has_value() );
        BOOST_TEST( !r.has_error() );

        BOOST_TEST_EQ( std::move( r ).error().v_, 0 );
    }

    {
        double x = 1.0;

        BOOST_TEST(( result<double&, X>( x ).has_value() ));
        BOOST_TEST(( !result<double&, X>( x ).has_error() ));

        BOOST_TEST_EQ( (result<double&, X>( x ).error().v_), 0 );
    }

    //

    return boost::report_errors();
}
