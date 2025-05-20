// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#if defined(__GNUC__) && __GNUC__ == 7
# pragma GCC diagnostic ignored "-Wnoexcept-type"
#endif

#include <boost/compat/invoke.hpp>
#include <boost/config.hpp>
#include <boost/config/workaround.hpp>

#define STATIC_ASSERT(...) static_assert(__VA_ARGS__, #__VA_ARGS__)
#define BOOST_TEST_EQ(x, y) STATIC_ASSERT((x) == (y))

struct X
{
    constexpr int f0() const
    {
        return -1;
    }

    constexpr int f1( int x1 ) const noexcept
    {
        return x1;
    }

    constexpr int f2( int x1, int x2 ) const
    {
        return 10*x1+x2;
    }

    constexpr int f3( int x1, int x2, int x3 ) const noexcept
    {
        return 100*x1 + 10*x2 + x3;
    }
};

struct Y: public X
{
};

int main()
{
    {
        BOOST_TEST_EQ( boost::compat::invoke( &X::f0, X() ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::f1, X(), 1 ), 1 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::f2, X(), 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::f3, X(), 1, 2, 3 ), 123 );
    }

    {
        BOOST_TEST_EQ( boost::compat::invoke( &X::f0, Y() ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::f1, Y(), 1 ), 1 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::f2, Y(), 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::f3, Y(), 1, 2, 3 ), 123 );
    }

    {
        constexpr X x = {};

        BOOST_TEST_EQ( boost::compat::invoke( &X::f0, x ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::f1, x, 1 ), 1 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::f2, x, 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::f3, x, 1, 2, 3 ), 123 );

        BOOST_TEST_EQ( boost::compat::invoke( &X::f0, &x ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::f1, &x, 1 ), 1 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::f2, &x, 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::f3, &x, 1, 2, 3 ), 123 );
    }

#if !BOOST_WORKAROUND(BOOST_MSVC, < 1910)
    {
        constexpr Y y = {};

        BOOST_TEST_EQ( boost::compat::invoke( &X::f0, y ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::f1, y, 1 ), 1 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::f2, y, 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::f3, y, 1, 2, 3 ), 123 );

        BOOST_TEST_EQ( boost::compat::invoke( &X::f0, &y ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::f1, &y, 1 ), 1 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::f2, &y, 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::f3, &y, 1, 2, 3 ), 123 );
    }
#endif
}
