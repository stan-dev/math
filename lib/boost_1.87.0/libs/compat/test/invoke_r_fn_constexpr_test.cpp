// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#if defined(__GNUC__) && __GNUC__ == 7
# pragma GCC diagnostic ignored "-Wnoexcept-type"
#endif

#include <boost/compat/invoke.hpp>
#include <boost/config.hpp>

#define STATIC_ASSERT(...) static_assert((__VA_ARGS__), #__VA_ARGS__)
#define BOOST_TEST_EQ(x, y) STATIC_ASSERT((x) == (y))

constexpr int f0()
{
    return -1;
}

constexpr int f1( int x1 ) noexcept
{
    return x1;
}

constexpr int f2( int x1, int x2 )
{
    return 10*x1+x2;
}

constexpr int f3( int x1, int x2, int x3 ) noexcept
{
    return 100*x1 + 10*x2 + x3;
}

int main()
{
    BOOST_TEST_EQ( boost::compat::invoke_r<long>( f0 ), -1 );
    BOOST_TEST_EQ( boost::compat::invoke_r<long>( f1, 1 ), 1 );
    BOOST_TEST_EQ( boost::compat::invoke_r<long>( f2, 1, 2 ), 12 );
    BOOST_TEST_EQ( boost::compat::invoke_r<long>( f3, 1, 2, 3 ), 123 );

#if !defined(BOOST_NO_CXX14_CONSTEXPR)

    STATIC_ASSERT( boost::compat::invoke_r<void>( f0 ), true );
    STATIC_ASSERT( boost::compat::invoke_r<void>( f1, 1 ), true );
    STATIC_ASSERT( boost::compat::invoke_r<void>( f2, 1, 2 ), true );
    STATIC_ASSERT( boost::compat::invoke_r<void>( f3, 1, 2, 3 ), true );

#endif
}
