// Copyright 2017, 2018, 2024 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/detail/sha1.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/core/snprintf.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cstdint>
#include <string>
#include <cstddef>

template<std::size_t N> std::string digest_to_string( unsigned char const (&v)[ N ] )
{
    std::string r;

    for( std::size_t i = 0; i < N; ++i )
    {
        char buffer[ 8 ];
        boost::core::snprintf( buffer, sizeof( buffer ), "%02x", static_cast<int>( v[ i ] ) );

        r += buffer;
    }

    return r;
}

template<class Hash> std::string digest( std::string const& s )
{
    Hash hash;
    hash.process_bytes( s.data(), s.size() );

    typename Hash::digest_type result;
    hash.get_digest( result );

    return digest_to_string( result );
}

int main()
{
    using boost::uuids::detail::sha1;

    // Test vectors from https://en.wikipedia.org/wiki/SHA-1

    BOOST_TEST_EQ( digest<sha1>( "" ), std::string( "da39a3ee5e6b4b0d3255bfef95601890afd80709" ) );
    BOOST_TEST_EQ( digest<sha1>( "The quick brown fox jumps over the lazy dog" ), std::string( "2fd4e1c67a2d28fced849ee1bb76e7391b93eb12" ) );
    BOOST_TEST_EQ( digest<sha1>( "The quick brown fox jumps over the lazy cog" ), std::string( "de9f2c7fd25e1b3afad3e85a0bd17d9b100db4b3" ) );

    // Test vectors from https://tools.ietf.org/html/rfc3174

    BOOST_TEST_EQ( digest<sha1>( "abc" ), std::string( "a9993e364706816aba3e25717850c26c9cd0d89d" ) );
    BOOST_TEST_EQ( digest<sha1>( "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq" ), std::string( "84983e441c3bd26ebaae4aa1f95129e5e54670f1" ) );
    BOOST_TEST_EQ( digest<sha1>( std::string( 1000000, 'a' ) ), std::string( "34aa973cd4c4daa4f61eeb2bdbad27316534016f" ) );

    using boost::uuids::uuid;

    BOOST_TEST_EQ( sha1().get_version(), uuid::version_name_based_sha1 );

    return boost::report_errors();
}
