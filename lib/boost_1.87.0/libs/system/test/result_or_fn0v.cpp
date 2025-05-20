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

int f()
{
    return 2;
}

X g()
{
    return { 2 };
}

int& h()
{
    static int x = 2;
    return x;
}

int main()
{
    {
        result<int> r( 1 );

        int x = r | f;

        BOOST_TEST_EQ( x, 1 );
    }

    {
        result<int> const r( 1 );

        int x = r | f;

        BOOST_TEST_EQ( x, 1 );
    }

    {
        int x = result<int>( 1 ) | f;

        BOOST_TEST_EQ( x, 1 );
    }

    {
        result<int, E> r( in_place_error );

        int x = r | f;

        BOOST_TEST_EQ( x, 2 );
    }

    {
        result<int, E> const r( in_place_error );

        int x = r | f;

        BOOST_TEST_EQ( x, 2 );
    }

    {
        int x = result<int, E>( in_place_error ) | f;

        BOOST_TEST_EQ( x, 2 );
    }

    {
        Y y = result<Y>( in_place_value, 1 ) | g;

        BOOST_TEST_EQ( y.v_, 1 );
    }

    {
        Y y = result<Y, E>( in_place_error ) | g;

        BOOST_TEST_EQ( y.v_, 2 );
    }

    {
        int x1 = 1;

        result<int&> r( x1 );

        int& x = r | h;

        BOOST_TEST_EQ( &x, &x1 );
    }

    {
        int x1 = 1;

        result<int&> const r( x1 );

        int& x = r | h;

        BOOST_TEST_EQ( &x, &x1 );
    }

    {
        int x1 = 1;

        int& x = result<int&>( x1 ) | h;

        BOOST_TEST_EQ( &x, &x1 );
    }

    {
        result<int&, E> r( in_place_error );

        int& x = r | h;

        BOOST_TEST_EQ( &x, &h() );
    }

    {
        result<int&, E> const r( in_place_error );

        int& x = r | h;

        BOOST_TEST_EQ( &x, &h() );
    }

    {
        int& x = result<int&, E>( in_place_error ) | h;

        BOOST_TEST_EQ( &x, &h() );
    }

    return boost::report_errors();
}
