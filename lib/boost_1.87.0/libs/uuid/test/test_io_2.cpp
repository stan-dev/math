// Copyright 2009 Andy Tompkins
// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/core/lightweight_test.hpp>
#include <string>
#include <sstream>
#include <iomanip>

using namespace boost::uuids;

std::string narrow( std::string const& in )
{
    return in;
}

std::string narrow( std::wstring const& in )
{
    std::string out( in.size(), 0 );

    for( std::size_t i = 0; i < in.size(); ++i )
    {
        out[ i ] = static_cast<char>( static_cast<unsigned char>( in[ i ] ) );
    }

    return out;
}

template<class Ch> void uufail( std::basic_string<Ch> const& in )
{
    std::basic_istringstream<Ch> is( in );

    uuid out;

    if( !BOOST_TEST_NOT( is >> out ) )
    {
        std::cerr << "  in='" << narrow( in ) << "'\n";
    }
}

template<class Ch, class T1> void uusucceed( Ch const* in, T1 const& t1 )
{
    std::basic_istringstream<Ch> is( in );

    {
        T1 out;
        is >> out;

        BOOST_TEST_EQ( is.fail(), false );
        BOOST_TEST_EQ( out, t1 );
    }
}

template<class Ch, class T1, class T2> void uusucceed( Ch const* in, T1 const& t1, T2 const& t2 )
{
    std::basic_istringstream<Ch> is( in );

    {
        T1 out;
        is >> out;

        BOOST_TEST_EQ( is.fail(), false );
        BOOST_TEST_EQ( out, t1 );
    }

    {
        T2 out;
        is >> out;

        BOOST_TEST_EQ( is.fail(), false );
        BOOST_TEST_EQ( out, t2 );
    }
}

template<class Ch, class T1, class T2, class T3> void uusucceed( Ch const* in, T1 const& t1, T2 const& t2, T3 const& t3 )
{
    std::basic_istringstream<Ch> is( in );

    {
        T1 out;
        is >> out;

        BOOST_TEST_EQ( is.fail(), false );
        BOOST_TEST_EQ( out, t1 );
    }

    {
        T2 out;
        is >> out;

        BOOST_TEST_EQ( is.fail(), false );
        BOOST_TEST_EQ( out, t2 );
    }

    {
        T3 out;
        is >> out;

        BOOST_TEST_EQ( is.fail(), false );
        BOOST_TEST_EQ( out, t3 );
    }
}

int main()
{
    const uuid u1 = {{0}};
    const uuid u2 = {{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}};
    const uuid u3 = {{0x12,0x34,0x56,0x78,0x90,0xab,0xcd,0xef,0x12,0x34,0x56,0x78,0x90,0xab,0xcd,0xef}};

    uusucceed( "00000000-0000-0000-0000-000000000000", u1 );
    uusucceed( "00010203-0405-0607-0809-0a0b0c0d0e0f", u2 );
    uusucceed( "12345678-90ab-cdef-1234-567890abcdef", u3 );

    uusucceed( "00010203-0405-0607-0809-0A0B0C0D0E0F", u2 );
    uusucceed( "12345678-90AB-CDEF-1234-567890ABCDEF", u3 );

    uusucceed( "12345678-90AB-CDEF-1234-567890abcdef", u3 );
    uusucceed( "12345678-90ab-cdef-1234-567890ABCDEF", u3 );

    uusucceed( L"00000000-0000-0000-0000-000000000000", u1 );
    uusucceed( L"00010203-0405-0607-0809-0a0b0c0d0e0f", u2 );
    uusucceed( L"12345678-90ab-cdef-1234-567890abcdef", u3 );

    uusucceed( L"00010203-0405-0607-0809-0A0B0C0D0E0F", u2 );
    uusucceed( L"12345678-90AB-CDEF-1234-567890ABCDEF", u3 );

    uusucceed( L"12345678-90AB-CDEF-1234-567890abcdef", u3 );
    uusucceed( L"12345678-90ab-cdef-1234-567890ABCDEF", u3 );

    uusucceed( " 00000000-0000-0000-0000-000000000000", u1 );
    uusucceed( "  00010203-0405-0607-0809-0a0b0c0d0e0f", u2 );
    uusucceed( "   12345678-90AB-CDEF-1234-567890ABCDEF", u3 );

    uusucceed( L" 00000000-0000-0000-0000-000000000000", u1 );
    uusucceed( L"  00010203-0405-0607-0809-0a0b0c0d0e0f", u2 );
    uusucceed( L"   12345678-90AB-CDEF-1234-567890ABCDEF", u3 );

    uusucceed( " 00000000-0000-0000-0000-000000000000 ", u1 );
    uusucceed( "  00010203-0405-0607-0809-0A0B0C0D0E0F  ", u2 );
    uusucceed( "   12345678-90ab-cdef-1234-567890abcdef   ", u3 );

    uusucceed( L" 00000000-0000-0000-0000-000000000000 ", u1 );
    uusucceed( L"  00010203-0405-0607-0809-0A0B0C0D0E0F  ", u2 );
    uusucceed( L"   12345678-90ab-cdef-1234-567890abcdef   ", u3 );

    uusucceed( "00000000-0000-0000-0000-000000000000 00010203-0405-0607-0809-0a0b0c0d0e0f", u1, u2 );
    uusucceed( " 00000000-0000-0000-0000-000000000000  00010203-0405-0607-0809-0A0B0C0D0E0F", u1, u2 );
    uusucceed( "  00000000-0000-0000-0000-000000000000   00010203-0405-0607-0809-0a0b0c0d0e0f   ", u1, u2 );

    uusucceed( L"00000000-0000-0000-0000-000000000000 00010203-0405-0607-0809-0a0b0c0d0e0f", u1, u2 );
    uusucceed( L" 00000000-0000-0000-0000-000000000000  00010203-0405-0607-0809-0A0B0C0D0E0F", u1, u2 );
    uusucceed( L"  00000000-0000-0000-0000-000000000000   00010203-0405-0607-0809-0a0b0c0d0e0f   ", u1, u2 );

    uusucceed( "00000000-0000-0000-0000-000000000000 00010203-0405-0607-0809-0a0b0c0d0e0f 12345678-90AB-CDEF-1234-567890ABCDEF", u1, u2, u3 );
    uusucceed( " 00000000-0000-0000-0000-000000000000  00010203-0405-0607-0809-0A0B0C0D0E0F   12345678-90ab-cdef-1234-567890abcdef", u1, u2, u3 );
    uusucceed( "  00000000-0000-0000-0000-000000000000   00010203-0405-0607-0809-0a0b0c0d0e0f    12345678-90AB-CDEF-1234-567890ABCDEF     ", u1, u2, u3 );

    uusucceed( L"00000000-0000-0000-0000-000000000000 00010203-0405-0607-0809-0A0B0C0D0E0F 12345678-90ab-cdef-1234-567890abcdef", u1, u2, u3 );
    uusucceed( L" 00000000-0000-0000-0000-000000000000  00010203-0405-0607-0809-0a0b0c0d0e0f   12345678-90AB-CDEF-1234-567890ABCDEF", u1, u2, u3 );
    uusucceed( L"  00000000-0000-0000-0000-000000000000   00010203-0405-0607-0809-0A0B0C0D0E0F    12345678-90ab-cdef-1234-567890abcdef     ", u1, u2, u3 );

    uusucceed( "00010203-0405-0607-0809-0a0b0c0d0e0f -11 12345678-90AB-CDEF-1234-567890ABCDEF", u2, -11, u3 );
    uusucceed( L"00010203-0405-0607-0809-0a0b0c0d0e0f -11 12345678-90AB-CDEF-1234-567890ABCDEF", u2, -11, u3 );

    uusucceed( " 00010203-0405-0607-0809-0A0B0C0D0E0F  A    12345678-90ab-cdef-1234-567890abcdef    ", u2, 'A', u3 );
    uusucceed( L" 00010203-0405-0607-0809-0A0B0C0D0E0F  A    12345678-90ab-cdef-1234-567890abcdef    ", u2, L'A', u3 );

    std::string s2( "00010203-0405-0607-0809-0a0b0c0d0e0f" );
    std::wstring w2( L"00010203-0405-0607-0809-0a0b0c0d0e0f" );

    for( std::size_t i = 0; i < s2.size(); ++i )
    {
        uufail( s2.substr( 0, i ) );
    }

    for( std::size_t i = 0; i < w2.size(); ++i )
    {
        uufail( w2.substr( 0, i ) );
    }

    for( std::size_t i = 0; i < s2.size(); ++i )
    {
        std::string q2( s2 );

        q2[ i ] = 'g';
        uufail( q2 );
    }

    for( std::size_t i = 0; i < w2.size(); ++i )
    {
        std::wstring q2( w2 );

        q2[ i ] = L'g';
        uufail( q2 );
    }

    for( std::size_t i = 0; i < s2.size(); ++i )
    {
        std::string q2( s2 );

        if( q2[ i ] == '-' )
        {
            q2[ i ] = 'a';
        }
        else
        {
            q2[ i ] = '-';
        }

        uufail( q2 );
    }

    for( std::size_t i = 0; i < w2.size(); ++i )
    {
        std::wstring q2( w2 );

        if( q2[ i ] == L'-' )
        {
            q2[ i ] = L'a';
        }
        else
        {
            q2[ i ] = L'-';
        }

        uufail( q2 );
    }

    return boost::report_errors();
}
