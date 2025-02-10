// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/name_generator_sha1.hpp>
#include <boost/uuid/string_generator.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::uuids;

void test( uuid const& nmsp, char const* name, char const* expected )
{
    name_generator_sha1 gen( nmsp );

    uuid u0 = string_generator()( expected );

    uuid u1 = gen( name );
    BOOST_TEST_EQ( u1, u0 );

    uuid u2 = gen( std::string( name ) );
    BOOST_TEST_EQ( u2, u0 );

    uuid u3 = gen( name, std::strlen( name ) );
    BOOST_TEST_EQ( u3, u0 );
}

template<class Ch> void test( uuid const& nmsp, Ch const* name, char const* expected )
{
    name_generator_sha1 gen( nmsp );

    uuid u0 = string_generator()( expected );

    uuid u1 = gen( name );
    BOOST_TEST_EQ( u1, u0 );

    uuid u2 = gen( std::basic_string<Ch>( name ) );
    BOOST_TEST_EQ( u2, u0 );
}

int main()
{
    // https://stackoverflow.com/questions/5515880/test-vectors-for-uuid-version-5-converting-hash-into-guid-generation-algorith
    // https://de.wikipedia.org/wiki/Universally_Unique_Identifier#Namensbasierte_UUIDs_.28Version_3_und_5.29

    test( ns::dns(), "www.example.org", "74738ff5-5367-5958-9aee-98fffdcd1876" );
    test( ns::dns(), u"www.example.org", "74738ff5-5367-5958-9aee-98fffdcd1876" );
    test( ns::dns(), U"www.example.org", "74738ff5-5367-5958-9aee-98fffdcd1876" );
    test( ns::dns(), u8"www.example.org", "74738ff5-5367-5958-9aee-98fffdcd1876" );

    // RFC 4122bis Section A.4

    test( ns::dns(), "www.example.com", "2ed6657d-e927-568b-95e1-2665a8aea6a2" );
    test( ns::dns(), u"www.example.com", "2ed6657d-e927-568b-95e1-2665a8aea6a2" );
    test( ns::dns(), U"www.example.com", "2ed6657d-e927-568b-95e1-2665a8aea6a2" );
    test( ns::dns(), u8"www.example.com", "2ed6657d-e927-568b-95e1-2665a8aea6a2" );

    // https://uuid.ramsey.dev/en/stable/rfc4122/version5.html

    test( ns::url(), "https://www.php.net", "a8f6ae40-d8a7-58f0-be05-a22f94eca9ec" );
    test( ns::url(), u"https://www.php.net", "a8f6ae40-d8a7-58f0-be05-a22f94eca9ec" );
    test( ns::url(), U"https://www.php.net", "a8f6ae40-d8a7-58f0-be05-a22f94eca9ec" );
    test( ns::url(), u8"https://www.php.net", "a8f6ae40-d8a7-58f0-be05-a22f94eca9ec" );

    // test cases from test_name_generator.cpp

    test( ns::dns(), "www.widgets.com", "21f7f8de-8051-5b89-8680-0195ef798b6a" );
    test( ns::dns(), L"www.widgets.com", "c315270b-a466-5872-aca4-9626cec0f4be" );

    // examples from documentation

    test( ns::dns(), "boost.org", "0043f363-bbb4-5369-840a-322df6ec1926" );
    test( ns::dns(), L"boost.org", "c31c5016-3493-5dc2-8484-5813d495cc18" );

    // test wide strings

    {
        uuid nmsp = string_generator()( "70a4abc5-80ab-4176-8e11-bc5836b6fef9" );

        name_generator_sha1 gen( nmsp );

        wchar_t const name[] = { 1, 2, 3, 4, 0 };
        unsigned char const name2[] = { 1, 0, 0, 0, 2, 0, 0, 0, 3, 0, 0, 0, 4, 0, 0, 0 };

        uuid u1 = gen( name );
        uuid u2 = gen( name2, sizeof( name2 ) );

        BOOST_TEST_EQ( u1, u2 );
    }

    // test unicode strings

    {
        uuid nmsp = string_generator()( "70a4abc5-80ab-4176-8e11-bc5836b6fef9" );

        name_generator_sha1 gen( nmsp );

        char32_t const name32[] = { 0x0024, 0x00A3, 0x0418, 0x0939, 0x20AC, 0xD55C, 0xDC12, 0xD834, 0x10348, 0x1096B3, 0xD956, 0 };
        char16_t const name16[] = { 0x0024, 0x00A3, 0x0418, 0x0939, 0x20AC, 0xD55C, 0xDC12, 0xD834, 0xD800, 0xDF48, 0xDBE5, 0xDEB3, 0xD956, 0 };

        unsigned char name8[] =
        {
            /*U+0024*/ 0x24,
            /*U+00A3*/ 0xC2, 0xA3,
            /*U+0418*/ 0xD0, 0x98,
            /*U+0939*/ 0xE0, 0xA4, 0xB9,
            /*U+20AC*/ 0xE2, 0x82, 0xAC,
            /*U+D55C*/ 0xED, 0x95, 0x9C,
            /*U+DC12*/ 0xEF, 0xBF, 0xBD, // U+FFFD
            /*U+D834*/ 0xEF, 0xBF, 0xBD, // U+FFFD
            /*U+10348*/ 0xF0, 0x90, 0x8D, 0x88,
            /*U+1096B3*/ 0xF4, 0x89, 0x9A, 0xB3,
            /*U+D956*/ 0xEF, 0xBF, 0xBD, // U+FFFD
        };

        uuid u1 = gen( name32 );
        uuid u2 = gen( name16 );
        uuid u3 = gen( name8, sizeof( name8 ) );

        BOOST_TEST_EQ( u1, u3 );
        BOOST_TEST_EQ( u2, u3 );
    }

    return boost::report_errors();
}
