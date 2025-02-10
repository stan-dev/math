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
    BOOST_TEST_EQ( noexcept( boost::compat::invoke( f0 ) ), false );
    BOOST_TEST_EQ( noexcept( boost::compat::invoke( f1, 1 ) ), true );
    BOOST_TEST_EQ( noexcept( boost::compat::invoke( f2, 1, 2 ) ), false );
    BOOST_TEST_EQ( noexcept( boost::compat::invoke( f3, 1, 2, 3 ) ), true );

    return boost::report_errors();
}

#endif
