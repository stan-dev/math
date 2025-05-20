// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#if defined(__GNUC__) && __GNUC__ == 7
# pragma GCC diagnostic ignored "-Wnoexcept-type"
#endif

#include <boost/compat/bind_front.hpp>
#include <boost/core/lightweight_test.hpp>
#include <memory>

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

int g( std::unique_ptr<int> p, std::unique_ptr<int> q )
{
    return 10 * *p + *q;
}

int main()
{
    BOOST_TEST_EQ( boost::compat::bind_front( f0 )(), -1 );
    BOOST_TEST_EQ( boost::compat::bind_front( f1 )( 1 ), 1 );
    BOOST_TEST_EQ( boost::compat::bind_front( f2 )( 1, 2 ), 12 );
    BOOST_TEST_EQ( boost::compat::bind_front( f3 )( 1, 2, 3 ), 123 );

    BOOST_TEST_EQ( boost::compat::bind_front( f1, 1 )(), 1 );
    BOOST_TEST_EQ( boost::compat::bind_front( f2, 1 )( 2 ), 12 );
    BOOST_TEST_EQ( boost::compat::bind_front( f3, 1 )( 2, 3 ), 123 );

    BOOST_TEST_EQ( boost::compat::bind_front( f2, 1, 2 )(), 12 );
    BOOST_TEST_EQ( boost::compat::bind_front( f3, 1, 2 )( 3 ), 123 );

    BOOST_TEST_EQ( boost::compat::bind_front( f3, 1, 2, 3 )(), 123 );

    BOOST_TEST_EQ( boost::compat::bind_front( g, std::unique_ptr<int>( new int(1) ) )( std::unique_ptr<int>( new int(2) ) ), 12 );

    return boost::report_errors();
}
