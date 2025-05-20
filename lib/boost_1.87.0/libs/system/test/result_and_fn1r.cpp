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

struct E2
{
    E2() {}
    E2( E ) {}
};

result<int, E2> fi( int x )
{
    return 2 * x + 1;
}

result<int, E2> fi2( int )
{
    return E2();
}

result<void, E2> fi3( int )
{
    return {};
}

result<X, E2> fy( Y y )
{
    return X{ 2 * y.v_ + 1 };
}

result<X, E2> fy2( Y )
{
    return E2();
}

result<void, E2> fy3( Y )
{
    return {};
}

result<int, E2> fri( int& x )
{
    return x * 2 + 1;
}

result<int&, E2> fri2( int& )
{
    return E2();
}

result<void, E2> fri3( int& )
{
    return {};
}

result<int, E2> fk()
{
    return 7;
}

result<int, E2> fk2()
{
    return E2();
}

result<void, E2> fk3()
{
    return {};
}

result<void, E2> fk4()
{
    return E2();
}

int main()
{
    {
        result<int, E> r( 1 );

        {
            result<int, E2> r2 = r & fi;
            BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( *r2, 3 );
        }

        {
            result<int, E2> r2 = r & fi2;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<void, E2> r2 = r & fi3;
            BOOST_TEST( r2.has_value() );
        }
    }

    {
        result<int, E> const r( 1 );

        {
            result<int, E2> r2 = r & fi;
            BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( *r2, 3 );
        }

        {
            result<int, E2> r2 = r & fi2;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<void, E2> r2 = r & fi3;
            BOOST_TEST( r2.has_value() );
        }
    }

    {
        result<int, E2> r2 = result<int, E>( 1 ) & fi;
        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( *r2, 3 );
    }

    {
        result<int, E2> r2 = result<int, E>( 1 ) & fi2;
        BOOST_TEST( r2.has_error() );
    }

    {
        result<void, E2> r2 = result<int, E>( 1 ) & fi3;
        BOOST_TEST( r2.has_value() );
    }

    {
        result<int, E> r( in_place_error );

        {
            result<int, E2> r2 = r & fi;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<int, E2> r2 = r & fi2;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<void, E2> r2 = r & fi3;
            BOOST_TEST( r2.has_error() );
        }
    }

    {
        result<int, E> const r( in_place_error );

        {
            result<int, E2> r2 = r & fi;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<int, E2> r2 = r & fi2;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<void, E2> r2 = r & fi3;
            BOOST_TEST( r2.has_error() );
        }
    }

    {
        result<int, E2> r2 = result<int, E>( in_place_error ) & fi;
        BOOST_TEST( r2.has_error() );
    }

    {
        result<int, E2> r2 = result<int, E>( in_place_error ) & fi2;
        BOOST_TEST( r2.has_error() );
    }

    {
        result<void, E2> r2 = result<int, E>( in_place_error ) & fi3;
        BOOST_TEST( r2.has_error() );
    }

    {
        result<X, E2> r2 = result<Y, E>( in_place_value, 1 ) & fy;

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( r2->v_, 3 );
    }

    {
        result<X, E2> r2 = result<Y, E>( in_place_value, 1 ) & fy2;

        BOOST_TEST( r2.has_error() );
    }

    {
        result<void, E2> r2 = result<Y, E>( in_place_value, 1 ) & fy3;

        BOOST_TEST( r2.has_value() );
    }

    {
        result<X, E2> r2 = result<Y, E>( in_place_error ) & fy;

        BOOST_TEST( r2.has_error() );
    }

    {
        result<X, E2> r2 = result<Y, E>( in_place_error ) & fy2;

        BOOST_TEST( r2.has_error() );
    }

    {
        result<void, E2> r2 = result<Y, E>( in_place_error ) & fy3;

        BOOST_TEST( r2.has_error() );
    }

    {
        int x1 = 1;
        result<int&, E> r( x1 );

        {
            result<int, E2> r2 = r & fri;
            BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( *r2, 3 );
        }

        {
            result<int&, E2> r2 = r & fri2;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<void, E2> r2 = r & fri3;
            BOOST_TEST( r2.has_value() );
        }
    }

    {
        int x1 = 1;
        result<int&, E> const r( x1 );

        {
            result<int, E2> r2 = r & fri;
            BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( *r2, 3 );
        }

        {
            result<int&, E2> r2 = r & fri2;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<void, E2> r2 = r & fri3;
            BOOST_TEST( r2.has_value() );
        }
    }

    {
        int x1 = 1;

        result<int, E2> r2 = result<int&, E>( x1 ) & fri;
        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( *r2, 3 );
    }

    {
        int x1 = 1;

        result<int&, E2> r2 = result<int&, E>( x1 ) & fri2;
        BOOST_TEST( r2.has_error() );
    }

    {
        int x1 = 1;

        result<void, E2> r2 = result<int&, E>( x1 ) & fri3;
        BOOST_TEST( r2.has_value() );
    }

    {
        result<int&, E> r( in_place_error );

        {
            result<int, E2> r2 = r & fri;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<int&, E2> r2 = r & fri2;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<void, E2> r2 = r & fri3;
            BOOST_TEST( r2.has_error() );
        }
    }

    {
        result<int&, E> const r( in_place_error );

        {
            result<int, E2> r2 = r & fri;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<int&, E2> r2 = r & fri2;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<void, E2> r2 = r & fri3;
            BOOST_TEST( r2.has_error() );
        }
    }

    {
        result<int, E2> r2 = result<int&, E>( in_place_error ) & fri;
        BOOST_TEST( r2.has_error() );
    }

    {
        result<int&, E2> r2 = result<int&, E>( in_place_error ) & fri2;
        BOOST_TEST( r2.has_error() );
    }

    {
        result<void, E2> r2 = result<int&, E>( in_place_error ) & fri3;
        BOOST_TEST( r2.has_error() );
    }

    {
        result<void, E> r;

        {
            result<int, E2> r2 = r & fk;
            BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( *r2, 7 );
        }

        {
            result<int, E2> r2 = r & fk2;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<void, E2> r2 = r & fk3;
            BOOST_TEST( r2.has_value() );
        }

        {
            result<void, E2> r2 = r & fk4;
            BOOST_TEST( r2.has_error() );
        }
    }

    {
        result<void, E> const r;

        {
            result<int, E2> r2 = r & fk;
            BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( *r2, 7 );
        }

        {
            result<int, E2> r2 = r & fk2;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<void, E2> r2 = r & fk3;
            BOOST_TEST( r2.has_value() );
        }

        {
            result<void, E2> r2 = r & fk4;
            BOOST_TEST( r2.has_error() );
        }
    }

    {
        result<int, E2> r2 = result<void, E>() & fk;
        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( *r2, 7 );
    }

    {
        result<int, E2> r2 = result<void, E>() & fk2;
        BOOST_TEST( r2.has_error() );
    }

    {
        result<void, E2> r2 = result<void, E>() & fk3;
        BOOST_TEST( r2.has_value() );
    }

    {
        result<void, E2> r2 = result<void, E>() & fk4;
        BOOST_TEST( r2.has_error() );
    }

    {
        result<void, E> r( in_place_error );

        {
            result<int, E2> r2 = r & fk;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<int, E2> r2 = r & fk2;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<void, E2> r2 = r & fk3;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<void, E2> r2 = r & fk4;
            BOOST_TEST( r2.has_error() );
        }
    }

    {
        result<void, E> const r( in_place_error );

        {
            result<int, E2> r2 = r & fk;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<int, E2> r2 = r & fk2;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<void, E2> r2 = r & fk3;
            BOOST_TEST( r2.has_error() );
        }

        {
            result<void, E2> r2 = r & fk4;
            BOOST_TEST( r2.has_error() );
        }
    }

    {
        result<int, E2> r2 = result<void, E>( in_place_error ) & fk;
        BOOST_TEST( r2.has_error() );
    }

    {
        result<int, E2> r2 = result<void, E>( in_place_error ) & fk2;
        BOOST_TEST( r2.has_error() );
    }

    {
        result<void, E2> r2 = result<void, E>( in_place_error ) & fk3;
        BOOST_TEST( r2.has_error() );
    }

    {
        result<void, E2> r2 = result<void, E>( in_place_error ) & fk4;
        BOOST_TEST( r2.has_error() );
    }

    return boost::report_errors();
}
