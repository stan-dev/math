// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/random_generator.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/config.hpp>
#include <boost/config/pragma_message.hpp>
#include <random>
#include <vector>
#include <cstring>

#if defined(BOOST_UUID_HAS_THREE_WAY_COMPARISON)

BOOST_PRAGMA_MESSAGE( "Three way comparisons are tested: BOOST_UUID_HAS_THREE_WAY_COMPARISON=" BOOST_STRINGIZE(BOOST_UUID_HAS_THREE_WAY_COMPARISON) )

#elif !defined(__cpp_impl_three_way_comparison)

BOOST_PRAGMA_MESSAGE( "Three way comparisons not tested because __cpp_impl_three_way_comparison is not defined" )

#elif !defined(__cpp_lib_three_way_comparison)

BOOST_PRAGMA_MESSAGE( "Three way comparisons not tested because __cpp_lib_three_way_comparison is not defined" )

#else

BOOST_PRAGMA_MESSAGE( "Three way comparisons not tested: __cpp_impl_three_way_comparison=" BOOST_STRINGIZE(__cpp_impl_three_way_comparison) ", __cpp_lib_three_way_comparison=" BOOST_STRINGIZE(__cpp_lib_three_way_comparison) )

#endif

using namespace boost::uuids;

void test_comparison( uuid const& u1, uuid const& u2 )
{
    int r = std::memcmp( u1.data, u2.data, 16 );

    if( r == 0 )
    {
        BOOST_TEST_EQ( u1, u2 );
    }

    if( r < 0 )
    {
        BOOST_TEST_LT( u1, u2 );
    }

    if( r > 0 )
    {
        BOOST_TEST_GT( u1, u2 );
    }

    if( r <= 0 )
    {
        BOOST_TEST_LE( u1, u2 );
    }

    if( r >= 0 )
    {
        BOOST_TEST_GE( u1, u2 );
    }

#if defined(BOOST_UUID_HAS_THREE_WAY_COMPARISON)

    constexpr auto eq = std::strong_ordering::equal;
    constexpr auto lt = std::strong_ordering::less;
    constexpr auto gt = std::strong_ordering::greater;

    if( r == 0 )
    {
        BOOST_TEST( ( u1 <=> u2 ) == eq );
    }

    if( r < 0 )
    {
        BOOST_TEST( ( u1 <=> u2 ) == lt );
    }

    if( r > 0 )
    {
        BOOST_TEST( ( u1 <=> u2 ) == gt );
    }

#endif
}

int main()
{
    std::vector<uuid> v;

    {
        uuid u1 = {};
        v.push_back( u1 );

        for( int i = 0; i < 16; ++i )
        {
            uuid u2( u1 );
            u2.data[ i ] = 0x01;

            v.push_back( u2 );
        }
    }

    {
        uuid u1 = {{ 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF }};
        v.push_back( u1 );

        for( int i = 0; i < 16; ++i )
        {
            uuid u2( u1 );
            u2.data[ i ] = 0xFE;

            v.push_back( u2 );
        }
    }

    {
        std::mt19937 rng;
        basic_random_generator<std::mt19937> gen( &rng );

        for( int i = 0; i < 16; ++i )
        {
            v.push_back( gen() );
        }
    }

    for( std::size_t i = 0, n = v.size(); i < n; ++i )
    {
        for( std::size_t j = 0; j < n; ++j )
        {
            test_comparison( v[ i ], v[ j ] );
        }
    }

    return boost::report_errors();
}
