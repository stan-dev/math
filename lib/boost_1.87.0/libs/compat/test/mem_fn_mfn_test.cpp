// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#if defined(__GNUC__) && __GNUC__ == 7
# pragma GCC diagnostic ignored "-Wnoexcept-type"
#endif

#include <boost/compat/mem_fn.hpp>
#include <boost/core/lightweight_test.hpp>
#include <functional>

struct X
{
    int f0()
    {
        return -1;
    }

    int f1( int x1 ) noexcept
    {
        return x1;
    }

    int f2( int x1, int x2 ) const
    {
        return 10*x1+x2;
    }

    int f3( int x1, int x2, int x3 ) const noexcept
    {
        return 100*x1 + 10*x2 + x3;
    }
};

struct Y: public virtual X
{
};

int main()
{
    {
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f0 )( X() ), -1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f1 )( X(), 1 ), 1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f2 )( X(), 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f3 )( X(), 1, 2, 3 ), 123 );
    }

    {
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f0 )( Y() ), -1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f1 )( Y(), 1 ), 1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f2 )( Y(), 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f3 )( Y(), 1, 2, 3 ), 123 );
    }

    {
        X x;

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f0 )( x ), -1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f1 )( x, 1 ), 1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f2 )( x, 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f3 )( x, 1, 2, 3 ), 123 );

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f0 )( &x ), -1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f1 )( &x, 1 ), 1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f2 )( &x, 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f3 )( &x, 1, 2, 3 ), 123 );

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f0 )( std::ref(x) ), -1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f1 )( std::ref(x), 1 ), 1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f2 )( std::ref(x), 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f3 )( std::ref(x), 1, 2, 3 ), 123 );

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f2 )( std::cref(x), 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f3 )( std::cref(x), 1, 2, 3 ), 123 );
    }

    {
        Y y;

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f0 )( y ), -1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f1 )( y, 1 ), 1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f2 )( y, 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f3 )( y, 1, 2, 3 ), 123 );

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f0 )( &y ), -1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f1 )( &y, 1 ), 1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f2 )( &y, 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f3 )( &y, 1, 2, 3 ), 123 );

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f0 )( std::ref(y) ), -1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f1 )( std::ref(y), 1 ), 1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f2 )( std::ref(y), 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f3 )( std::ref(y), 1, 2, 3 ), 123 );

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f2 )( std::cref(y), 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f3 )( std::cref(y), 1, 2, 3 ), 123 );
    }

    {
        X const x = {};

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f2 )( x, 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f3 )( x, 1, 2, 3 ), 123 );

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f2 )( &x, 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f3 )( &x, 1, 2, 3 ), 123 );

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f2 )( std::ref(x), 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f3 )( std::ref(x), 1, 2, 3 ), 123 );
    }

    {
        Y const y = {};

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f2 )( y, 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f3 )( y, 1, 2, 3 ), 123 );

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f2 )( &y, 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f3 )( &y, 1, 2, 3 ), 123 );

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f2 )( std::ref(y), 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::f3 )( std::ref(y), 1, 2, 3 ), 123 );
    }

    return boost::report_errors();
}
