// Copyright 2017, 2021, 2023 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/system/result.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/core/lightweight_test_trait.hpp>

using namespace boost::system;

struct X
{
    static int instances;

    int v_;

    explicit X( int v ): v_( v ) { ++instances; }

    X( X const& ) = delete;
    X& operator=( X const& ) = delete;

    ~X() { --instances; }
};

int X::instances = 0;

int main()
{
    {
        int x = 0;
        result<int&> r( x );

        BOOST_TEST( r.has_value() );
        BOOST_TEST( !r.has_error() );

        BOOST_TEST_EQ( r.value(), 0 );
    }

    {
        int x = 0;
        result<int&> r = x;

        BOOST_TEST( r.has_value() );
        BOOST_TEST( !r.has_error() );

        BOOST_TEST_EQ( r.value(), 0 );
    }

    {
        int x = 1;
        result<int&, int> r( in_place_value, x );

        BOOST_TEST( r.has_value() );
        BOOST_TEST( !r.has_error() );

        BOOST_TEST_EQ( *r, 1 );
    }

    BOOST_TEST_EQ( X::instances, 0 );

    {
        X x( 1 );
        result<X&> r( x );

        BOOST_TEST( r.has_value() );
        BOOST_TEST( !r.has_error() );

        BOOST_TEST_EQ( r.value().v_, 1 );

        BOOST_TEST_EQ( X::instances, 1 );
    }

    BOOST_TEST_EQ( X::instances, 0 );

    {
        X x( 1 );
        result<X&> r = x;

        BOOST_TEST( r.has_value() );
        BOOST_TEST( !r.has_error() );

        BOOST_TEST_EQ( r.value().v_, 1 );

        BOOST_TEST_EQ( X::instances, 1 );
    }

    BOOST_TEST_EQ( X::instances, 0 );

    {
        X x( 1 );
        result<X&, X> r( in_place_value, x );

        BOOST_TEST( r.has_value() );
        BOOST_TEST( !r.has_error() );

        BOOST_TEST_EQ( r->v_, 1 );

        BOOST_TEST_EQ( X::instances, 1 );
    }

    BOOST_TEST_EQ( X::instances, 0 );

    {
        BOOST_TEST_TRAIT_TRUE((std::is_constructible<result<int&>, int&>));
        BOOST_TEST_TRAIT_TRUE((std::is_convertible<int&, result<int&>>));

        BOOST_TEST_TRAIT_FALSE((std::is_constructible<result<int&>, int const&>));
        BOOST_TEST_TRAIT_FALSE((std::is_convertible<int const&, result<int&>>));

        BOOST_TEST_TRAIT_FALSE((std::is_constructible<result<int&>, int>));
        BOOST_TEST_TRAIT_FALSE((std::is_convertible<int, result<int&>>));

        BOOST_TEST_TRAIT_FALSE((std::is_constructible<result<int&, int>, int&>));
        BOOST_TEST_TRAIT_FALSE((std::is_constructible<result<int&, float>, int&>));

        BOOST_TEST_TRAIT_FALSE((std::is_constructible<result<X&>, int>));
        BOOST_TEST_TRAIT_FALSE((std::is_convertible<int, result<X&>>));
    }

    return boost::report_errors();
}
