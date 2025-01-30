// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/compat/mem_fn.hpp>
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
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( X() ), -1 );
    }

    {
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( Y() ), -1 );
    }

    {
        X x;

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( x ), -1 );
        boost::compat::mem_fn( &X::m )( x ) = +1;
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( x ), +1 );
    }

    {
        X x;

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( &x ), -1 );
        boost::compat::mem_fn( &X::m )( &x ) = +1;
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( &x ), +1 );
    }

    {
        X x;

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( std::ref(x) ), -1 );
        boost::compat::mem_fn( &X::m )( std::ref(x) ) = +1;
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( std::ref(x) ), +1 );
    }

    {
        X x;

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( std::cref(x) ), -1 );
    }

    {
        Y y;

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( y ), -1 );
        boost::compat::mem_fn( &X::m )( y ) = +1;
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( y ), +1 );
    }

    {
        Y y;

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( &y ), -1 );
        boost::compat::mem_fn( &X::m )( &y ) = +1;
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( &y ), +1 );
    }

    {
        Y y;

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( std::ref(y) ), -1 );
        boost::compat::mem_fn( &X::m )( std::ref(y) ) = +1;
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( std::ref(y) ), +1 );
    }

    {
        Y y;

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( std::cref(y) ), -1 );
    }

    {
        X const x = {};

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( x ), -1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( &x ), -1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( std::ref(x) ), -1 );
    }

    {
        Y const y = {};

        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( y ), -1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( &y ), -1 );
        BOOST_TEST_EQ( boost::compat::mem_fn( &X::m )( std::ref(y) ), -1 );
    }

    return boost::report_errors();
}
