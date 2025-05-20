// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/compat/bind_back.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/config.hpp>
#include <boost/config/workaround.hpp>

struct F1
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

struct F2
{
    int operator()( int x1, int x2 ) &
    {
        return 100*x1 + 10*x2 + 1;
    }

    int operator()( int x1, int x2 ) const &
    {
        return 100*x1 + 10*x2 + 2;
    }
    int operator()( int x1, int x2 ) &&
    {
        return 100*x1 + 10*x2 + 3;
    }
    int operator()( int x1, int x2 ) const &&
    {
        return 100*x1 + 10*x2 + 4;
    }
};

int main()
{
    {
        BOOST_TEST_EQ( boost::compat::bind_back( F1() )(), -1 );
        BOOST_TEST_EQ( boost::compat::bind_back( F1() )( 1 ), 1 );
        BOOST_TEST_EQ( boost::compat::bind_back( F1() )( 1, 2 ), 12 );
        BOOST_TEST_EQ( boost::compat::bind_back( F1() )( 1, 2, 3 ), 123 );

        BOOST_TEST_EQ( boost::compat::bind_back( F1(), 1 )(), 1 );
        BOOST_TEST_EQ( boost::compat::bind_back( F1(), 1 )( 2 ), 21 );
        BOOST_TEST_EQ( boost::compat::bind_back( F1(), 1 )( 2, 3 ), 231 );

        BOOST_TEST_EQ( boost::compat::bind_back( F1(), 1, 2 )(), 12 );
        BOOST_TEST_EQ( boost::compat::bind_back( F1(), 1, 2 )( 3 ), 312 );

        BOOST_TEST_EQ( boost::compat::bind_back( F1(), 1, 2, 3 )(), 123 );
    }

    {
        BOOST_TEST_EQ( boost::compat::bind_back( F2(), 9 )( 8 ), 893 );
    }

    {
        auto fn = boost::compat::bind_back( F2(), 9 );
        BOOST_TEST_EQ( fn( 8 ), 891 );
    }

    {
        auto const fn = boost::compat::bind_back( F2(), 9 );
        BOOST_TEST_EQ( fn( 8 ), 892 );
    }

    {
        auto fn = boost::compat::bind_back( F2(), 9 );
        BOOST_TEST_EQ( std::move( fn )( 8 ), 893 );
    }

#if !BOOST_WORKAROUND(BOOST_GCC, < 40900)

    {
        auto const fn = boost::compat::bind_back( F2(), 9 );
        BOOST_TEST_EQ( std::move( fn )( 8 ), 894 );
    }

#endif

    return boost::report_errors();
}
