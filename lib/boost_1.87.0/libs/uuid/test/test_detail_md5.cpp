// Copyright 2017, 2018, 2024 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/detail/md5.hpp>
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

template<class Hash> std::string digest( char const * s )
{
    Hash hash;
    hash.process_bytes( s, std::strlen( s ) );

    typename Hash::digest_type result;
    hash.get_digest( result );

    return digest_to_string( result );
}

int main()
{
    using boost::uuids::detail::md5;

    // Test vectors from https://tools.ietf.org/html/rfc1321

    BOOST_TEST_EQ( digest<md5>( "" ), std::string( "d41d8cd98f00b204e9800998ecf8427e" ) );
    BOOST_TEST_EQ( digest<md5>( "a" ), std::string( "0cc175b9c0f1b6a831c399e269772661" ) );
    BOOST_TEST_EQ( digest<md5>( "abc" ), std::string( "900150983cd24fb0d6963f7d28e17f72" ) );
    BOOST_TEST_EQ( digest<md5>( "message digest" ), std::string( "f96b697d7cb7938d525a2f31aaf161d0" ) );
    BOOST_TEST_EQ( digest<md5>( "abcdefghijklmnopqrstuvwxyz" ), std::string( "c3fcd3d76192e4007dfb496cca67e13b" ) );
    BOOST_TEST_EQ( digest<md5>( "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789" ), std::string( "d174ab98d277d9f5a5611c2c9f419d9f" ) );
    BOOST_TEST_EQ( digest<md5>( "12345678901234567890123456789012345678901234567890123456789012345678901234567890" ), std::string( "57edf4a22be3c955ac49da2e2107b67a" ) );

    using boost::uuids::uuid;

    BOOST_TEST_EQ( md5().get_version(), uuid::version_name_based_md5 );

    return boost::report_errors();
}
