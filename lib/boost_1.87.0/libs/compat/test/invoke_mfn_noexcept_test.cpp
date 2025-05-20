// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/config/pragma_message.hpp>

#if !defined(__cpp_noexcept_function_type)

BOOST_PRAGMA_MESSAGE( "Test skipped, __cpp_noexcept_function_type is not defined" )
int main() {}

#else

#include <boost/compat/invoke.hpp>
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
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f0, X() ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f1, X(), 1 ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f2, X(), 1, 2 ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f3, X(), 1, 2, 3 ) ), true );
    }

    {
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f0, Y() ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f1, Y(), 1 ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f2, Y(), 1, 2 ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f3, Y(), 1, 2, 3 ) ), true );
    }

    {
        X x;

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f0, x ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f1, x, 1 ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f2, x, 1, 2 ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f3, x, 1, 2, 3 ) ), true );

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f0, &x ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f1, &x, 1 ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f2, &x, 1, 2 ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f3, &x, 1, 2, 3 ) ), true );

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f0, std::ref(x) ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f1, std::ref(x), 1 ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f2, std::ref(x), 1, 2 ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f3, std::ref(x), 1, 2, 3 ) ), true );

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f2, std::cref(x), 1, 2 ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f3, std::cref(x), 1, 2, 3 ) ), true );
    }

    {
        Y y;

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f0, y ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f1, y, 1 ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f2, y, 1, 2 ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f3, y, 1, 2, 3 ) ), true );

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f0, &y ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f1, &y, 1 ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f2, &y, 1, 2 ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f3, &y, 1, 2, 3 ) ), true );

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f0, std::ref(y) ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f1, std::ref(y), 1 ) ), true );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f2, std::ref(y), 1, 2 ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f3, std::ref(y), 1, 2, 3 ) ), true );

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f2, std::cref(y), 1, 2 ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f3, std::cref(y), 1, 2, 3 ) ), true );
    }

    {
        X const x = {};

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f2, x, 1, 2 ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f3, x, 1, 2, 3 ) ), true );

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f2, &x, 1, 2 ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f3, &x, 1, 2, 3 ) ), true );

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f2, std::ref(x), 1, 2 ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f3, std::ref(x), 1, 2, 3 ) ), true );
    }

    {
        Y const y = {};

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f2, y, 1, 2 ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f3, y, 1, 2, 3 ) ), true );

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f2, &y, 1, 2 ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f3, &y, 1, 2, 3 ) ), true );

        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f2, std::ref(y), 1, 2 ) ), false );
        BOOST_TEST_EQ( noexcept( boost::compat::invoke( &X::f3, std::ref(y), 1, 2, 3 ) ), true );
    }

    return boost::report_errors();
}

#endif
