// Copyright 2017, 2021-2024 Peter Dimov.
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

result<int, E> fi()
{
    return 2;
}

result<int, E2> fi2()
{
    return E2();
}

result<X, E> fy()
{
    return X{ 2 };
}

result<Y, E2> fy2()
{
    return E2();
}

result<int&, E> fri()
{
    static int x = 2;
    return x;
}

result<int&, E2> fri2()
{
    return E2();
}

int main()
{
    {
        result<int, E2> r( 1 );

        r |= fi;

        BOOST_TEST( r.has_value() ) && BOOST_TEST_EQ( *r, 1 );

        r |= fi2;

        BOOST_TEST( r.has_value() ) && BOOST_TEST_EQ( *r, 1 );
    }

    {
        result<int, E2> r( in_place_error );

        r |= fi2;

        BOOST_TEST( r.has_error() );

        r |= fi;

        BOOST_TEST( r.has_value() ) && BOOST_TEST_EQ( *r, 2 );
    }

    {
        result<Y, E2> r( in_place_value, 1 );

        r |= fy;

        BOOST_TEST( r.has_value() ) && BOOST_TEST_EQ( r->v_, 1 );

        r |= fy2;

        BOOST_TEST( r.has_value() ) && BOOST_TEST_EQ( r->v_, 1 );
    }

    {
        result<Y, E2> r( in_place_error );

        r |= fy2;

        BOOST_TEST( r.has_error() );

        r |= fy;

        BOOST_TEST( r.has_value() ) && BOOST_TEST_EQ( r->v_, 2 );
    }

    {
        int x1 = 1;

        result<int&, E2> r( x1 );

        r |= fri;

        BOOST_TEST( r.has_value() ) && BOOST_TEST_EQ( &*r, &x1 );

        r |= fri2;

        BOOST_TEST( r.has_value() ) && BOOST_TEST_EQ( &*r, &x1 );
    }

    {
        result<int&, E2> r( in_place_error );

        r |= fri2;

        BOOST_TEST( r.has_error() );

        r |= fri;

        BOOST_TEST( r.has_value() ) && BOOST_TEST_EQ( &*r, &*fri() );
    }

    return boost::report_errors();
}
