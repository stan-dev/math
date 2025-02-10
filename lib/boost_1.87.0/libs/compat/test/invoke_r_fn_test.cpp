// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#if defined(__GNUC__) && __GNUC__ == 7
# pragma GCC diagnostic ignored "-Wnoexcept-type"
#endif

#include <boost/compat/invoke.hpp>
#include <boost/core/lightweight_test.hpp>

int f0()
{
    return -1;
}

int f1( int x1 ) noexcept
{
    return x1;
}

int f2( int x1, int x2 )
{
    return 10*x1+x2;
}

int f3( int x1, int x2, int x3 ) noexcept
{
    return 100*x1 + 10*x2 + x3;
}

int main()
{
    BOOST_TEST_EQ( boost::compat::invoke_r<long>( f0 ), -1 );
    BOOST_TEST_EQ( boost::compat::invoke_r<long>( f1, 1 ), 1 );
    BOOST_TEST_EQ( boost::compat::invoke_r<long>( f2, 1, 2 ), 12 );
    BOOST_TEST_EQ( boost::compat::invoke_r<long>( f3, 1, 2, 3 ), 123 );

    boost::compat::invoke_r<void>( f0 );
    boost::compat::invoke_r<void>( f1, 1 );
    boost::compat::invoke_r<void>( f2, 1, 2 );
    boost::compat::invoke_r<void>( f3, 1, 2, 3 );

    return boost::report_errors();
}
