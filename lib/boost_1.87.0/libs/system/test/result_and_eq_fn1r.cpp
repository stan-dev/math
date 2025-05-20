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

    Y& operator=( Y&& r )
    {
        if( &r != this )
        {
            v_ = r.v_;
            r.v_ = 0;
        }

        return *this;
    }
};

struct E
{
};

struct E2
{
    E2() {}
    E2( E ) {}
};

result<int, E> fi( int x )
{
    return 2 * x + 1;
}

result<int, E2> fi2( int )
{
    return E2();
}

result<X, E> fy( Y y )
{
    return X{ 2 * y.v_ + 1 };
}

result<Y, E2> fy2( Y )
{
    return E2();
}

result<int&, E> fri( int& )
{
    static int x = 2;
    return x;
}

result<int&, E2> fri2( int& )
{
    return E2();
}

result<void, E2> fv()
{
    return {};
}

result<void, E2> fv2()
{
    return E2();
}

int main()
{
    {
        result<int, E2> r( 1 );

        r &= fi;
        BOOST_TEST( r.has_value() ) && BOOST_TEST_EQ( *r, 3 );

        r &= fi2;
        BOOST_TEST( r.has_error() );

        r &= fi;
        BOOST_TEST( r.has_error() );
    }

    {
        result<Y, E2> r( in_place_value, 1 );

        r &= fy;
        BOOST_TEST( r.has_value() ) && BOOST_TEST_EQ( r->v_, 3 );

        r &= fy2;
        BOOST_TEST( r.has_error() );

        r &= fy;
        BOOST_TEST( r.has_error() );
    }

    {
        int x1 = 1;
        result<int&, E2> r( x1 );

        r &= fri;
        BOOST_TEST( r.has_value() ) && BOOST_TEST_EQ( &*r, &*fri( x1 ) );

        r &= fri2;
        BOOST_TEST( r.has_error() );

        r &= fri;
        BOOST_TEST( r.has_error() );
    }

    {
        result<void, E2> r;

        r &= fv;
        BOOST_TEST( r.has_value() );

        r &= fv2;
        BOOST_TEST( r.has_error() );

        r &= fv;
        BOOST_TEST( r.has_error() );
    }

    return boost::report_errors();
}
