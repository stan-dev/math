// Copyright 2017, 2021, 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/system/result.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::system;

struct X
{
    static int instances;

    int v_;

    explicit X( int v ): v_( v ) { ++instances; }

    X( int v1, int v2 ): v_( v1+v2 ) { ++instances; }
    X( int v1, int v2, int v3 ): v_( v1+v2+v3 ) { ++instances; }

    X( X const& ) = delete;
    X& operator=( X const& ) = delete;

    ~X() { --instances; }
};

int X::instances = 0;

int main()
{
    {
        using R = result<int>;
        R r( R::in_place_value, 0 );

        BOOST_TEST( r.has_value() );
        BOOST_TEST( !r.has_error() );

        BOOST_TEST_EQ( r.value(), 0 );
    }

    {
        using R = result<int, int>;
        R r( R::in_place_value, 1 );

        BOOST_TEST( r.has_value() );
        BOOST_TEST( !r.has_error() );

        BOOST_TEST_EQ( *r, 1 );
    }

    BOOST_TEST_EQ( X::instances, 0 );

    {
        using R = result<X>;
        R r( R::in_place_value, 1 );

        BOOST_TEST( r.has_value() );
        BOOST_TEST( !r.has_error() );

        BOOST_TEST_EQ( r.value().v_, 1 );

        BOOST_TEST_EQ( X::instances, 1 );
    }

    BOOST_TEST_EQ( X::instances, 0 );

    {
        using R = result<X>;
        R r( R::in_place_value, 1, 2 );

        BOOST_TEST( r.has_value() );
        BOOST_TEST( !r.has_error() );

        BOOST_TEST_EQ( r.value().v_, 1+2 );

        BOOST_TEST_EQ( X::instances, 1 );
    }

    BOOST_TEST_EQ( X::instances, 0 );

    {
        using R = result<X>;
        R r( R::in_place_value, 1, 2, 3 );

        BOOST_TEST( r.has_value() );
        BOOST_TEST( !r.has_error() );

        BOOST_TEST_EQ( r.value().v_, 1+2+3 );

        BOOST_TEST_EQ( X::instances, 1 );
    }

    BOOST_TEST_EQ( X::instances, 0 );

    {
        using R = result<X, X>;
        R r( R::in_place_value, 1 );

        BOOST_TEST( r.has_value() );
        BOOST_TEST( !r.has_error() );

        BOOST_TEST_EQ( r->v_, 1 );

        BOOST_TEST_EQ( X::instances, 1 );
    }

    BOOST_TEST_EQ( X::instances, 0 );

    {
        using R = result<void>;
        R r( R::in_place_value );

        BOOST_TEST( r.has_value() );
        BOOST_TEST( !r.has_error() );
    }

    return boost::report_errors();
}
