// Copyright 2017, 2021, 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/system/result.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cerrno>

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

struct Y
{
    static int instances;

    Y() noexcept { ++instances; }

    Y( Y const& ) noexcept { ++instances; }

    Y& operator=( Y const& ) = default;

    ~Y() { --instances; }
};

int Y::instances = 0;

int main()
{
    {
        result<int> r;

        BOOST_TEST( r.has_value() );

        r.emplace( 1 );

        BOOST_TEST( r.has_value() );
        BOOST_TEST_EQ( r.value(), 1 );
    }

    {
        result<int> r( ENOENT, generic_category() );

        BOOST_TEST( !r.has_value() );

        r.emplace( 1 );

        BOOST_TEST( r.has_value() );
        BOOST_TEST_EQ( r.value(), 1 );
    }

    BOOST_TEST_EQ( X::instances, 0 );

    {
        result<X> r( 0 );

        BOOST_TEST( r.has_value() );
        BOOST_TEST_EQ( r.value().v_, 0 );
        BOOST_TEST_EQ( X::instances, 1 );

        r.emplace( 1 );
        BOOST_TEST( r.has_value() );
        BOOST_TEST_EQ( r.value().v_, 1 );
        BOOST_TEST_EQ( X::instances, 1 );

        r.emplace( 1, 2 );
        BOOST_TEST( r.has_value() );
        BOOST_TEST_EQ( r.value().v_, 1+2 );
        BOOST_TEST_EQ( X::instances, 1 );

        r.emplace( 1, 2, 3 );
        BOOST_TEST( r.has_value() );
        BOOST_TEST_EQ( r.value().v_, 1+2+3 );
        BOOST_TEST_EQ( X::instances, 1 );
    }

    BOOST_TEST_EQ( X::instances, 0 );

    {
        result<X> r( ENOENT, generic_category() );

        BOOST_TEST( !r.has_value() );
        BOOST_TEST_EQ( X::instances, 0 );

        r.emplace( 1, 2 );
        BOOST_TEST( r.has_value() );
        BOOST_TEST_EQ( r.value().v_, 1+2 );
        BOOST_TEST_EQ( X::instances, 1 );
    }

    BOOST_TEST_EQ( X::instances, 0 );
    BOOST_TEST_EQ( Y::instances, 0 );

    {
        result<int, Y> r( Y{} );

        BOOST_TEST( !r.has_value() );
        BOOST_TEST_EQ( Y::instances, 1 );

        r.emplace( 1 );
        BOOST_TEST( r.has_value() );
        BOOST_TEST_EQ( *r, 1 );
        BOOST_TEST_EQ( Y::instances, 0 );
    }

    BOOST_TEST_EQ( X::instances, 0 );
    BOOST_TEST_EQ( Y::instances, 0 );

    {
        result<X, Y> r( in_place_error );

        BOOST_TEST( !r.has_value() );
        BOOST_TEST_EQ( X::instances, 0 );
        BOOST_TEST_EQ( Y::instances, 1 );

        r.emplace( 1, 2, 3 );

        BOOST_TEST( r.has_value() );
        BOOST_TEST_EQ( r->v_, 1+2+3 );
        BOOST_TEST_EQ( X::instances, 1 );
        BOOST_TEST_EQ( Y::instances, 0 );
    }

    BOOST_TEST_EQ( X::instances, 0 );
    BOOST_TEST_EQ( Y::instances, 0 );

    {
        result<void> r;

        BOOST_TEST( r.has_value() );

        r.emplace();
        BOOST_TEST( r.has_value() );
    }

    {
        result<void> r( ENOENT, generic_category() );

        BOOST_TEST( !r.has_value() );

        r.emplace();
        BOOST_TEST( r.has_value() );
    }

    BOOST_TEST_EQ( Y::instances, 0 );

    {
        result<void, Y> r( in_place_error );

        BOOST_TEST( !r.has_value() );
        BOOST_TEST_EQ( Y::instances, 1 );

        r.emplace();
        BOOST_TEST( r.has_value() );
        BOOST_TEST_EQ( Y::instances, 0 );
    }

    return boost::report_errors();
}
