// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/compat/bind_back.hpp>
#include <boost/config.hpp>
#include <boost/config/pragma_message.hpp>

#if defined(BOOST_NO_CXX14_CONSTEXPR)

BOOST_PRAGMA_MESSAGE( "Test skipped because BOOST_NO_CXX14_CONSTEXPR is defined" )
int main() {}

#else

#define STATIC_ASSERT(...) static_assert(__VA_ARGS__, #__VA_ARGS__)
#define BOOST_TEST_EQ(x, y) STATIC_ASSERT((x) == (y))

struct F1
{
    constexpr int operator()()
    {
        return -1;
    }

    constexpr int operator()( int x1 ) noexcept
    {
        return x1;
    }

    constexpr int operator()( int x1, int x2 ) const
    {
        return 10*x1+x2;
    }

    constexpr int operator()( int x1, int x2, int x3 ) const noexcept
    {
        return 100*x1 + 10*x2 + x3;
    }
};

struct F2
{
    constexpr int operator()( int x1, int x2 ) &
    {
        return 100*x1 + 10*x2 + 1;
    }

    constexpr int operator()( int x1, int x2 ) const &
    {
        return 100*x1 + 10*x2 + 2;
    }
    constexpr int operator()( int x1, int x2 ) &&
    {
        return 100*x1 + 10*x2 + 3;
    }
    constexpr int operator()( int x1, int x2 ) const &&
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
        constexpr auto fn = boost::compat::bind_back( F2(), 9 );
        BOOST_TEST_EQ( fn( 8 ), 892 );
    }

    {
        constexpr auto fn = boost::compat::bind_back( F2(), 9 );
        BOOST_TEST_EQ( std::move( fn )( 8 ), 894 );
    }
}

#endif
