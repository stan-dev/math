// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/compat/invoke.hpp>
#include <boost/config.hpp>
#include <boost/config/workaround.hpp>

#define STATIC_ASSERT(...) static_assert(__VA_ARGS__, #__VA_ARGS__)
#define BOOST_TEST_EQ(x, y) STATIC_ASSERT((x) == (y))

struct X
{
    int m = -1;
};

struct Y: public X
{
};

int main()
{
    {
        BOOST_TEST_EQ( boost::compat::invoke( &X::m, X() ), -1 );

        constexpr X x = {};

        BOOST_TEST_EQ( boost::compat::invoke( &X::m, x ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::m, &x ), -1 );
    }

#if !BOOST_WORKAROUND(BOOST_MSVC, >= 1920 && BOOST_MSVC < 1950)
    {
        BOOST_TEST_EQ( boost::compat::invoke( &X::m, Y() ), -1 );

        constexpr Y y = {};

        BOOST_TEST_EQ( boost::compat::invoke( &X::m, y ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke( &X::m, &y ), -1 );
    }
#endif
}
