// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/compat/invoke.hpp>
#include <boost/config.hpp>
#include <boost/config/workaround.hpp>

#define STATIC_ASSERT(...) static_assert((__VA_ARGS__), #__VA_ARGS__)
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
        constexpr X x = {};

        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, X() ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, x ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, &x ), -1 );

#if !defined(BOOST_NO_CXX14_CONSTEXPR)

        STATIC_ASSERT( boost::compat::invoke_r<void>( &X::m, X() ), true );
        STATIC_ASSERT( boost::compat::invoke_r<void>( &X::m, x ), true );
        STATIC_ASSERT( boost::compat::invoke_r<void>( &X::m, &x ), true );

#endif
    }

#if !BOOST_WORKAROUND(BOOST_MSVC, >= 1920 && BOOST_MSVC < 1950)
    {
        constexpr Y y = {};

        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, Y() ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, y ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, &y ), -1 );

#if !defined(BOOST_NO_CXX14_CONSTEXPR)

        STATIC_ASSERT( boost::compat::invoke_r<void>( &X::m, Y() ), true );
        STATIC_ASSERT( boost::compat::invoke_r<void>( &X::m, y ), true );
        STATIC_ASSERT( boost::compat::invoke_r<void>( &X::m, &y ), true );

#endif
    }
#endif
}
