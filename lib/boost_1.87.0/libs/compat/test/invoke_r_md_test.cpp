// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/compat/invoke.hpp>
#include <boost/core/lightweight_test.hpp>
#include <functional>

struct X
{
    int m = -1;
};

struct Y: public virtual X
{
};

int main()
{
    {
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, X() ), -1 );
    }

    {
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, Y() ), -1 );
    }

    {
        X x;

        BOOST_TEST_EQ( boost::compat::invoke_r<int const&>( &X::m, x ), -1 );
        boost::compat::invoke_r<int&>( &X::m, x ) = +1;
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, x ), +1 );
        boost::compat::invoke_r<void>( &X::m, x );
    }

    {
        X x;

        BOOST_TEST_EQ( boost::compat::invoke_r<int const&>( &X::m, &x ), -1 );
        boost::compat::invoke_r<int&>( &X::m, &x ) = +1;
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, &x ), +1 );
        boost::compat::invoke_r<void>( &X::m, &x );
    }

    {
        X x;

        BOOST_TEST_EQ( boost::compat::invoke_r<int const&>( &X::m, std::ref(x) ), -1 );
        boost::compat::invoke_r<int&>( &X::m, std::ref(x) ) = +1;
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, std::ref(x) ), +1 );
        boost::compat::invoke_r<void>( &X::m, std::ref(x) );
    }

    {
        X x;

        BOOST_TEST_EQ( boost::compat::invoke_r<int const&>( &X::m, std::cref(x) ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, std::cref(x) ), -1 );
        boost::compat::invoke_r<void>( &X::m, std::cref(x) );
    }

    {
        Y y;

        BOOST_TEST_EQ( boost::compat::invoke_r<int const&>( &X::m, y ), -1 );
        boost::compat::invoke_r<int&>( &X::m, y ) = +1;
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, y ), +1 );
        boost::compat::invoke_r<void>( &X::m, y );
    }

    {
        Y y;

        BOOST_TEST_EQ( boost::compat::invoke_r<int const&>( &X::m, &y ), -1 );
        boost::compat::invoke_r<int&>( &X::m, &y ) = +1;
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, &y ), +1 );
        boost::compat::invoke_r<void>( &X::m, &y );
    }

    {
        Y y;

        BOOST_TEST_EQ( boost::compat::invoke_r<int const&>( &X::m, std::ref(y) ), -1 );
        boost::compat::invoke_r<int&>( &X::m, std::ref(y) ) = +1;
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, std::ref(y) ), +1 );
        boost::compat::invoke_r<void>( &X::m, std::ref(y) );
    }

    {
        Y y;

        BOOST_TEST_EQ( boost::compat::invoke_r<int const&>( &X::m, std::cref(y) ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, std::cref(y) ), -1 );
        boost::compat::invoke_r<void>( &X::m, std::cref(y) );
    }

    {
        X const x = {};

        BOOST_TEST_EQ( boost::compat::invoke_r<int const&>( &X::m, x ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, x ), -1 );
        boost::compat::invoke_r<void>( &X::m, x );

        BOOST_TEST_EQ( boost::compat::invoke_r<int const&>( &X::m, &x ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, &x ), -1 );
        boost::compat::invoke_r<void>( &X::m, &x );

        BOOST_TEST_EQ( boost::compat::invoke_r<int const&>( &X::m, std::ref(x) ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, std::ref(x) ), -1 );
        boost::compat::invoke_r<void>( &X::m, std::ref(x) );
    }

    {
        Y const y = {};

        BOOST_TEST_EQ( boost::compat::invoke_r<int const&>( &X::m, y ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, y ), -1 );
        boost::compat::invoke_r<void>( &X::m, y );

        BOOST_TEST_EQ( boost::compat::invoke_r<int const&>( &X::m, &y ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, &y ), -1 );
        boost::compat::invoke_r<void>( &X::m, &y );

        BOOST_TEST_EQ( boost::compat::invoke_r<int const&>( &X::m, std::ref(y) ), -1 );
        BOOST_TEST_EQ( boost::compat::invoke_r<long>( &X::m, std::ref(y) ), -1 );
        boost::compat::invoke_r<void>( &X::m, std::ref(y) );
    }

    return boost::report_errors();
}
