// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/charconv.hpp>
#include <boost/core/detail/splitmix64.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cstdio>

static boost::detail::splitmix64 rng;

template<class T> void zero_extend_test()
{
    int const N = 2048;
    int const M = 16;

    for( int i = 0; i < M; ++i )
    {
        char buffer[ static_cast<std::size_t>(N + 128) ];

        {
            unsigned long long v = static_cast<unsigned long long>( rng() );

            std::snprintf( buffer, sizeof( buffer ), "%llu", v );

            T ref;

            auto r1 = boost::charconv::from_chars( buffer, buffer + std::strlen( buffer ), ref );

            if( !BOOST_TEST_EQ( static_cast<int>(r1.ec), 0 ) )
            {
                std::fprintf( stderr, "Test failure for '%s': got error at position %ld\n", buffer, static_cast<long>( r1.ptr - buffer ) ); // LCOV_EXCL_LINE
            }
            else
            {
                for( int j = 1; j <= N; ++j )
                {
                    std::snprintf( buffer, sizeof( buffer ), "%llu%se-%d", v, std::string( static_cast<std::size_t>(j), '0' ).c_str(), j );

                    T w;
                    auto r2 = boost::charconv::from_chars( buffer, buffer + std::strlen( buffer ), w );

                    if( !BOOST_TEST_EQ( static_cast<int>(r2.ec), 0 ) )
                    {
                        std::fprintf( stderr, "Test failure for '%s': expected '%.15g', got error at position %ld\n", buffer, ref, static_cast<long>( r1.ptr - buffer ) ); // LCOV_EXCL_LINE
                    }
                    else if( !BOOST_TEST_EQ( w, ref ) )
                    {
                        std::fprintf( stderr, "Test failure for '%s': expected '%.15g', got '%.15g'\n", buffer, ref, w ); // LCOV_EXCL_LINE
                    }
                }

                for( int j = 1; j <= N; ++j )
                {
                    std::snprintf( buffer, sizeof( buffer ), "%llu0e-%s1", v, std::string( static_cast<std::size_t>(j), '0' ).c_str() );

                    T w;
                    auto r2 = boost::charconv::from_chars( buffer, buffer + std::strlen( buffer ), w );

                    if( !BOOST_TEST_EQ( static_cast<int>(r2.ec), 0 ) )
                    {
                        std::fprintf( stderr, "Test failure for '%s': expected '%.15g', got error at position %ld\n", buffer, ref, static_cast<long>( r1.ptr - buffer ) ); // LCOV_EXCL_LINE
                    }
                    else if( !BOOST_TEST_EQ( w, ref ) )
                    {
                        std::fprintf( stderr, "Test failure for '%s': expected '%.15g', got '%.15g'\n", buffer, ref, w ); // LCOV_EXCL_LINE
                    }
                }
            }
        }

        {
            unsigned long long v = static_cast<unsigned long long>( rng() );

            std::snprintf( buffer, sizeof( buffer ), "0.%llu", v );

            T ref;

            auto r1 = boost::charconv::from_chars( buffer, buffer + std::strlen( buffer ), ref );

            if( !BOOST_TEST_EQ( static_cast<int>(r1.ec), 0 ) )
            {
                std::fprintf( stderr, "Test failure for '%s': got error at position %ld\n", buffer, static_cast<long>( r1.ptr - buffer ) ); // LCOV_EXCL_LINE
            }
            else
            {
                for( int j = 1; j <= N; ++j )
                {
                    std::snprintf( buffer, sizeof( buffer ), "0.%s%llue%d", std::string( static_cast<std::size_t>(j), '0' ).c_str(), v, j );

                    T w;
                    auto r2 = boost::charconv::from_chars( buffer, buffer + std::strlen( buffer ), w );

                    if( !BOOST_TEST_EQ( static_cast<int>(r2.ec), 0 ) )
                    {
                        std::fprintf( stderr, "Test failure for '%s': expected '%.15g', got error at position %ld\n", buffer, ref, static_cast<long>( r1.ptr - buffer ) ); // LCOV_EXCL_LINE
                    }
                    else if( !BOOST_TEST_EQ( w, ref ) )
                    {
                        std::fprintf( stderr, "Test failure for '%s': expected '%.15g', got '%.15g'\n", buffer, ref, w ); // LCOV_EXCL_LINE
                    }
                }
            }
        }
    }
}

int main()
{
    zero_extend_test<float>();
    zero_extend_test<double>();

    return boost::report_errors();
}
