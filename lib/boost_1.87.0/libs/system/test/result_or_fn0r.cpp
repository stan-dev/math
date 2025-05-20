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

result<int, E> fi()
{
    return 2;
}

result<int, E> fi2()
{
    return E();
}

result<Y, E> fy()
{
    return Y{ 2 };
}

result<Y, E> fy2()
{
    return E();
}

result<int&, E> fri()
{
    static int x = 2;
    return x;
}

result<int&, E> fri2()
{
    return E();
}

result<void, E> fv()
{
    return {};
}

result<void, E> fv2()
{
    return E();
}

int main()
{
    {
        result<int> r( 1 );

        int x = r | fi | 3;

        BOOST_TEST_EQ( x, 1 );
    }

    {
        result<int> const r( 1 );

        int x = r | fi | 3;

        BOOST_TEST_EQ( x, 1 );
    }

    {
        int x = result<int>( 1 ) | fi | 3;

        BOOST_TEST_EQ( x, 1 );
    }

    {
        result<int> r( 1 );

        int x = r | fi2 | 3;

        BOOST_TEST_EQ( x, 1 );
    }

    {
        result<int> const r( 1 );

        int x = r | fi2 | 3;

        BOOST_TEST_EQ( x, 1 );
    }

    {
        int x = result<int>( 1 ) | fi2 | 3;

        BOOST_TEST_EQ( x, 1 );
    }

    {
        result<int> r( in_place_error );

        int x = r | fi | 3;

        BOOST_TEST_EQ( x, 2 );
    }

    {
        result<int> const r( in_place_error );

        int x = r | fi | 3;

        BOOST_TEST_EQ( x, 2 );
    }

    {
        int x = result<int>( in_place_error ) | fi | 3;

        BOOST_TEST_EQ( x, 2 );
    }

    {
        result<int> r( in_place_error );

        int x = r | fi2 | 3;

        BOOST_TEST_EQ( x, 3 );
    }

    {
        result<int> const r( in_place_error );

        int x = r | fi2 | 3;

        BOOST_TEST_EQ( x, 3 );
    }

    {
        int x = result<int>( in_place_error ) | fi2 | 3;

        BOOST_TEST_EQ( x, 3 );
    }

    {
        Y y = result<X>( X{1} ) | fy | X{3};

        BOOST_TEST_EQ( y.v_, 1 );
    }

    {
        Y y = result<X>( X{1} ) | fy2 | X{3};

        BOOST_TEST_EQ( y.v_, 1 );
    }

    {
        Y y = result<X, E>( in_place_error ) | fy | X{3};

        BOOST_TEST_EQ( y.v_, 2 );
    }

    {
        Y y = result<X, E>( in_place_error ) | fy2 | Y{3};

        BOOST_TEST_EQ( y.v_, 3 );
    }

    {
        int x1 = 1;
        int x3 = 3;

        result<int&> r( x1 );

        int& x = r | fri | x3;

        BOOST_TEST_EQ( &x, &x1 );
    }

    {
        int x1 = 1;
        int x3 = 3;

        result<int&> const r( x1 );

        int& x = r | fri | x3;

        BOOST_TEST_EQ( &x, &x1 );
    }

    {
        int x1 = 1;
        int x3 = 3;

        int& x = result<int&>( x1 ) | fri | x3;

        BOOST_TEST_EQ( &x, &x1 );
    }

    {
        int x1 = 1;
        int x3 = 3;

        result<int&> r( x1 );

        int& x = r | fri2 | x3;

        BOOST_TEST_EQ( &x, &x1 );
    }

    {
        int x1 = 1;
        int x3 = 3;

        result<int&> const r( x1 );

        int& x = r | fri2 | x3;

        BOOST_TEST_EQ( &x, &x1 );
    }

    {
        int x1 = 1;
        int x3 = 3;

        int& x = result<int&>( x1 ) | fri2 | x3;

        BOOST_TEST_EQ( &x, &x1 );
    }

    {
        int x3 = 3;

        result<int&, E> r( in_place_error );

        int& x = r | fri | x3;

        BOOST_TEST_EQ( &x, &*fri() );
    }

    {
        int x3 = 3;

        result<int&, E> const r( in_place_error );

        int& x = r | fri | x3;

        BOOST_TEST_EQ( &x, &*fri() );
    }

    {
        int x3 = 3;

        int& x = result<int&, E>( in_place_error ) | fri | x3;

        BOOST_TEST_EQ( &x, &*fri() );
    }

    {
        int x3 = 3;

        result<int&, E> r( in_place_error );

        int& x = r | fri2 | x3;

        BOOST_TEST_EQ( &x, &x3 );
    }

    {
        int x3 = 3;

        result<int&, E> const r( in_place_error );

        int& x = r | fri2 | x3;

        BOOST_TEST_EQ( &x, &x3 );
    }

    {
        int x3 = 3;

        int& x = result<int&, E>( in_place_error ) | fri2 | x3;

        BOOST_TEST_EQ( &x, &x3 );
    }

    {
        result<void> r;
        result<void, E> r2 = r | fv;

        BOOST_TEST( r2.has_value() );
    }

    {
        result<void> r;
        result<void, E> r2 = r | fv2;

        BOOST_TEST( r2.has_value() );
    }

    {
        result<void> r( in_place_error );
        result<void, E> r2 = r | fv;

        BOOST_TEST( r2.has_value() );
    }

    {
        result<void> r( in_place_error );
        result<void, E> r2 = r | fv2;

        BOOST_TEST( r2.has_error() );
    }

    {
        result<void, E> r2 = result<void>() | fv;
        BOOST_TEST( r2.has_value() );
    }

    {
        result<void, E> r2 = result<void>() | fv2;
        BOOST_TEST( r2.has_value() );
    }

    {
        result<void, E> r2 = result<void>( in_place_error ) | fv;
        BOOST_TEST( r2.has_value() );
    }

    {
        result<void, E> r2 = result<void>( in_place_error ) | fv2;
        BOOST_TEST( r2.has_error() );
    }

    return boost::report_errors();
}
