// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/compat/bind_back.hpp>
#include <boost/config.hpp>
#include <boost/config/workaround.hpp>
#include <boost/config/pragma_message.hpp>

#if defined(BOOST_NO_CXX14_CONSTEXPR)

BOOST_PRAGMA_MESSAGE( "Test skipped because BOOST_NO_CXX14_CONSTEXPR is defined" )
int main() {}

#else

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
        BOOST_TEST_EQ( boost::compat::bind_back( &X::m, X() )(), -1 );
    }

    {
        constexpr X x = {};

        BOOST_TEST_EQ( boost::compat::bind_back( &X::m, x )(), -1 );
        BOOST_TEST_EQ( boost::compat::bind_back( &X::m, &x )(), -1 );
    }

#if !BOOST_WORKAROUND(BOOST_MSVC, >= 1920 && BOOST_MSVC < 1950)

    {
        BOOST_TEST_EQ( boost::compat::bind_back( &X::m, Y() )(), -1 );
    }

    {
        constexpr Y y = {};

        BOOST_TEST_EQ( boost::compat::bind_back( &X::m, y )(), -1 );
        BOOST_TEST_EQ( boost::compat::bind_back( &X::m, &y )(), -1 );
    }

#endif
}

#endif
