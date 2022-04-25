// Copyright 2022 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/unordered/detail/implementation.hpp>
#include <boost/core/bit.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/cstdint.hpp>

template<class Policy, class SizeT> void test( SizeT x )
{
    if( x <= 4 )
    {
        BOOST_TEST_EQ( Policy::new_bucket_count( x ), 4u );
    }
    else
    {
        BOOST_TEST_EQ( Policy::new_bucket_count( x ), boost::core::bit_ceil( x ) );
    }

    BOOST_TEST_EQ( Policy::prev_bucket_count( x ), boost::core::bit_floor( x ) );
}

int main()
{
    {
        typedef boost::uint64_t SizeT;
        typedef boost::unordered::detail::mix64_policy<SizeT> policy;

        for( SizeT i = 1; i < 200; ++i )
        {
            test<policy>( i );
        }

        for( int i = 8; i < 64; ++i )
        {
            SizeT x = SizeT( 1 ) << i;

            test<policy>( x - 1 );
            test<policy>( x );
            test<policy>( x + 1 );
        }
    }

    {
        typedef boost::uint32_t SizeT;
        typedef boost::unordered::detail::mix32_policy<SizeT> policy;

        for( SizeT i = 1; i < 200; ++i )
        {
            test<policy>( i );
        }

        for( int i = 8; i < 32; ++i )
        {
            SizeT x = SizeT( 1 ) << i;

            test<policy>( x - 1 );
            test<policy>( x );
            test<policy>( x + 1 );
        }
    }

    return boost::report_errors();
}
