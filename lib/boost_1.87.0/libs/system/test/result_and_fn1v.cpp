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

int k()
{
    return 7;
}

static int fv1_called_with;

void fv1( int x )
{
    fv1_called_with = x;
}

static int fv2_called;

void fv2()
{
    ++fv2_called;
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

    {
        result<void> r;
        result<int> r2 = r & k;

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( *r2, 7 );
    }

    {
        result<void> const r;
        result<int> r2 = r & k;

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( *r2, 7 );
    }

    {
        result<int> r2 = result<void>() & k;

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( *r2, 7 );
    }

    {
        result<void, E> r( in_place_error );
        result<int, E> r2 = r & k;

        BOOST_TEST( r2.has_error() );
    }

    {
        result<void, E> const r( in_place_error );
        result<int, E> r2 = r & k;

        BOOST_TEST( r2.has_error() );
    }

    {
        result<int, E> r2 = result<void, E>( in_place_error ) & k;

        BOOST_TEST( r2.has_error() );
    }

    {
        fv1_called_with = 0;

        result<int> r( 1 );
        result<void> r2 = r & fv1;

        BOOST_TEST( r2.has_value() );
        BOOST_TEST_EQ( fv1_called_with, 1 );
    }

    {
        result<int> const r( 1 );
        result<void> r2 = r & fv1;

        BOOST_TEST( r2.has_value() );
        BOOST_TEST_EQ( fv1_called_with, 1 );
    }

    {
        fv1_called_with = 0;

        result<void> r2 = result<int>( 1 ) & fv1;

        BOOST_TEST( r2.has_value() );
        BOOST_TEST_EQ( fv1_called_with, 1 );
    }

    {
        fv1_called_with = 0;

        result<int, E> r( in_place_error );
        result<void, E> r2 = r & fv1;

        BOOST_TEST( r2.has_error() );
        BOOST_TEST_EQ( fv1_called_with, 0 );
    }

    {
        fv1_called_with = 0;

        result<int, E> const r( in_place_error );
        result<void, E> r2 = r & fv1;

        BOOST_TEST( r2.has_error() );
        BOOST_TEST_EQ( fv1_called_with, 0 );
    }

    {
        fv1_called_with = 0;

        result<void, E> r2 = result<int, E>( in_place_error ) & fv1;

        BOOST_TEST( r2.has_error() );
        BOOST_TEST_EQ( fv1_called_with, 0 );
    }

    {
        result<void> r;
        fv2_called = 0;

        result<void> r2 = r & fv2;

        BOOST_TEST( r2.has_value() );
        BOOST_TEST_EQ( fv2_called, 1 );
    }

    {
        result<void> const r;
        fv2_called = 0;

        result<void> r2 = r & fv2;

        BOOST_TEST( r2.has_value() );
        BOOST_TEST_EQ( fv2_called, 1 );
    }

    {
        fv2_called = 0;

        result<void> r2 = result<void>() & fv2;

        BOOST_TEST( r2.has_value() );
        BOOST_TEST_EQ( fv2_called, 1 );
    }

    {
        result<void, E> r( in_place_error );
        fv2_called = 0;

        result<void, E> r2 = r & fv2;

        BOOST_TEST( r2.has_error() );
        BOOST_TEST_EQ( fv2_called, 0 );
    }

    {
        result<void, E> const r( in_place_error );
        fv2_called = 0;

        result<void, E> r2 = r & fv2;

        BOOST_TEST( r2.has_error() );
        BOOST_TEST_EQ( fv2_called, 0 );
    }

    {
        fv2_called = 0;

        result<void, E> r2 = result<void, E>( in_place_error ) & fv2;

        BOOST_TEST( r2.has_error() );
        BOOST_TEST_EQ( fv2_called, 0 );
    }

    return boost::report_errors();
}
