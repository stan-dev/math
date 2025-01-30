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
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::m, X() ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::m, Y() ) ), true );
    }

    {
        X x;

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::m, x ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::m, &x ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::m, std::ref(x) ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::m, std::cref(x) ) ), true );
    }

    {
        Y y;

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::m, y ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::m, &y ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::m, std::ref(y) ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::m, std::cref(y) ) ), true );
    }

    {
        X const x = {};

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::m, x ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::m, &x ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::m, std::ref(x) ) ), true );
    }

    {
        Y const y = {};

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::m, y ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::m, &y ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::m, std::ref(y) ) ), true );
    }

    return boost::report_errors();
}
