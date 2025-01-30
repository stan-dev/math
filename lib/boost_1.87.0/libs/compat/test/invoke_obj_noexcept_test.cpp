// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/compat/invoke.hpp>
#include <boost/core/lightweight_test.hpp>

struct F
{
    int operator()()
    {
        return -1;
    }

    int operator()( int x1 ) noexcept
    {
        return x1;
    }

    int operator()( int x1, int x2 ) const
    {
        return 10*x1+x2;
    }

    int operator()( int x1, int x2, int x3 ) const noexcept
    {
        return 100*x1 + 10*x2 + x3;
    }
};

int main()
{
    {
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( F() ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( F(), 1 ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( F(), 1, 2 ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( F(), 1, 2, 3 ) ), true );
    }

    {
        F f = {};

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( f ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( f, 1 ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( f, 1, 2 ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( f, 1, 2, 3 ) ), true );
    }

    {
        F const f = {};

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( f, 1, 2 ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( f, 1, 2, 3 ) ), true );
    }

    return boost::report_errors();
}
