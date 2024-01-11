// Copyright 2017, 2021, 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/system/result.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::system;

struct X
{
    int v_;
};

struct Y
{
    int v_;

    explicit Y( int v ): v_( v ) {}
    Y( X x ): v_( x.v_) {}

    Y( Y const& ) = delete;
    Y& operator=( Y const& ) = delete;

    Y( Y&& r ): v_( r.v_ )
    {
        r.v_ = 0;
    }

    Y& operator=( Y&& ) = delete;
};

struct E
{
};

int f( int x )
{
    return x * 2 + 1;
}

X g( Y y )
{
    return X{ y.v_ * 2 + 1 };
}

int& h( int& )
{
    static int x = 2;
    return x;
}

int main()
{
    {
        result<int> r( 1 );
        result<int> r2 = r & f;

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( *r2, 3 );
    }

    {
        result<int> const r( 1 );
        result<int> r2 = r & f;

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( *r2, 3 );
    }

    {
        result<int> r2 = result<int>( 1 ) & f;

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( *r2, 3 );
    }

    {
        result<int, E> r( in_place_error );
        result<int, E> r2 = r & f;

        BOOST_TEST( r2.has_error() );
    }

    {
        result<int, E> const r( in_place_error );
        result<int, E> r2 = r & f;

        BOOST_TEST( r2.has_error() );
    }

    {
        result<int, E> r2 = result<int, E>( in_place_error ) & f;

        BOOST_TEST( r2.has_error() );
    }

    {
        result<X> r2 = result<Y>( in_place_value, 1 ) & g;

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( r2->v_, 3 );
    }

    {
        result<X, E> r2 = result<Y, E>( in_place_error ) & g;

        BOOST_TEST( r2.has_error() );
    }

    {
        int x1 = 1;

        result<int&> r( x1 );

        result<int&> r2 = r & h;

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( &*r2, &h( x1 ) );
    }

    {
        int x1 = 1;

        result<int&> const r( x1 );

        result<int&> r2 = r & h;

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( &*r2, &h( x1 ) );
    }

    {
        int x1 = 1;

        result<int&> r2 = result<int&>( x1 ) & h;

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( &*r2, &h( x1 ) );
    }

    {
        result<int&, E> r( in_place_error );

        result<int&, E> r2 = r & h;

        BOOST_TEST( r2.has_error() );
    }

    {
        result<int&, E> const r( in_place_error );

        result<int&, E> r2 = r & h;

        BOOST_TEST( r2.has_error() );
    }

    {
        result<int&, E> r2 = result<int&, E>( in_place_error ) & h;

        BOOST_TEST( r2.has_error() );
    }

    return boost::report_errors();
}
