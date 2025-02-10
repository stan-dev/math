// Copyright 2009 Andy Tompkins 2009
// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/core/lightweight_test.hpp>
#include <string>
#include <algorithm>

using namespace boost::uuids;

template<std::size_t N, class Ch> void test2( uuid const& u, Ch const* expected )
{
    Ch buffer[ N ];

    std::fill( buffer + 0, buffer + N, '@' );

    // always succeeds, invalid uses are detected at compile time
    BOOST_TEST( to_chars( u, buffer ) == buffer + 36 );

    // check that the correct, zero-terminated, string is produced
    BOOST_TEST( buffer[ 36 ] == 0 ) && BOOST_TEST( std::basic_string<Ch>( buffer ) == std::basic_string<Ch>( expected ) );

    // make sure nothing is written beyond buffer[ 36 ]
    BOOST_TEST( std::basic_string<Ch>( buffer + 37, buffer + N ) == std::basic_string<Ch>( N - 37, '@' ) );
}

template<class Ch> void test( uuid const& u, Ch const* expected )
{
    test2<37>( u, expected );
    test2<38>( u, expected );
    test2<39>( u, expected );
    test2<40>( u, expected );
}

int main()
{
    uuid const u1 = {{ 0 }};
    uuid const u2 = {{ 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f }};
    uuid const u3 = {{ 0x12, 0x34, 0x56, 0x78, 0x90, 0xab, 0xcd, 0xef, 0x12, 0x34, 0x56, 0x78, 0x90, 0xab, 0xcd, 0xef }};

    test( u1, "00000000-0000-0000-0000-000000000000" );
    test( u2, "00010203-0405-0607-0809-0a0b0c0d0e0f" );
    test( u3, "12345678-90ab-cdef-1234-567890abcdef" );

    test( u1, L"00000000-0000-0000-0000-000000000000" );
    test( u2, L"00010203-0405-0607-0809-0a0b0c0d0e0f" );
    test( u3, L"12345678-90ab-cdef-1234-567890abcdef" );

    test( u1, u"00000000-0000-0000-0000-000000000000" );
    test( u2, u"00010203-0405-0607-0809-0a0b0c0d0e0f" );
    test( u3, u"12345678-90ab-cdef-1234-567890abcdef" );

    test( u1, U"00000000-0000-0000-0000-000000000000" );
    test( u2, U"00010203-0405-0607-0809-0a0b0c0d0e0f" );
    test( u3, U"12345678-90ab-cdef-1234-567890abcdef" );

    test( u1, u8"00000000-0000-0000-0000-000000000000" );
    test( u2, u8"00010203-0405-0607-0809-0a0b0c0d0e0f" );
    test( u3, u8"12345678-90ab-cdef-1234-567890abcdef" );

    return boost::report_errors();
}
