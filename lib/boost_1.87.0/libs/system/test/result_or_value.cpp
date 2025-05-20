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

int main()
{
    {
        result<int> r( 1 );

        int x = r | 2;

        BOOST_TEST_EQ( x, 1 );
    }

    {
        result<int> const r( 1 );

        int x = r | 2;

        BOOST_TEST_EQ( x, 1 );
    }

    {
        int x = result<int>( 1 ) | 2;

        BOOST_TEST_EQ( x, 1 );
    }

    {
        result<int, E> r( in_place_error );

        int x = r | 2;

        BOOST_TEST_EQ( x, 2 );
    }

    {
        result<int, E> const r( in_place_error );

        int x = r | 2;

        BOOST_TEST_EQ( x, 2 );
    }

    {
        int x = result<int, E>( in_place_error ) | 2;

        BOOST_TEST_EQ( x, 2 );
    }

    {
        Y y = result<Y>( in_place_value, 1 ) | Y{2};

        BOOST_TEST_EQ( y.v_, 1 );
    }

    {
        Y y = result<Y, E>( in_place_error ) | Y{2};

        BOOST_TEST_EQ( y.v_, 2 );
    }

    {
        Y y = result<Y>( in_place_value, 1 ) | X{2};

        BOOST_TEST_EQ( y.v_, 1 );
    }

    {
        Y y = result<Y, E>( in_place_error ) | X{2};

        BOOST_TEST_EQ( y.v_, 2 );
    }

    {
        int x1 = 1;
        int x2 = 2;

        result<int&> r( x1 );

        int& x = r | x2;

        BOOST_TEST_EQ( &x, &x1 );
    }

    {
        int x1 = 1;
        int x2 = 2;

        result<int&> const r( x1 );

        int& x = r | x2;

        BOOST_TEST_EQ( &x, &x1 );
    }

    {
        int x1 = 1;
        int x2 = 2;

        int& x = result<int&>( x1 ) | x2;

        BOOST_TEST_EQ( &x, &x1 );
    }

    {
        int x2 = 2;

        result<int&, E> r( in_place_error );

        int& x = r | x2;

        BOOST_TEST_EQ( &x, &x2 );
    }

    {
        int x2 = 2;

        result<int&, E> const r( in_place_error );

        int& x = r | x2;

        BOOST_TEST_EQ( &x, &x2 );
    }

    {
        int x2 = 2;

        int& x = result<int&, E>( in_place_error ) | x2;

        BOOST_TEST_EQ( &x, &x2 );
    }

    return boost::report_errors();
}
