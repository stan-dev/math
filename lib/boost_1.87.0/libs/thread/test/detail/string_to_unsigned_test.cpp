// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/thread/detail/string_to_unsigned.hpp>
#include <boost/core/lightweight_test.hpp>

int main()
{
    using boost::thread_detail::string_to_unsigned;

    unsigned v;

    BOOST_TEST_NOT( string_to_unsigned( "", v ) ) && BOOST_TEST_EQ( v, 0 );
    BOOST_TEST_NOT( string_to_unsigned( " ", v ) ) && BOOST_TEST_EQ( v, 0 );
    BOOST_TEST_NOT( string_to_unsigned( "+1", v ) ) && BOOST_TEST_EQ( v, 0 );
    BOOST_TEST_NOT( string_to_unsigned( "-1", v ) ) && BOOST_TEST_EQ( v, 0 );
    BOOST_TEST_NOT( string_to_unsigned( "abc", v ) ) && BOOST_TEST_EQ( v, 0 );

    BOOST_TEST( string_to_unsigned( "0", v ) ) && BOOST_TEST_EQ( v, 0 );
    BOOST_TEST( string_to_unsigned( "1", v ) ) && BOOST_TEST_EQ( v, 1 );
    BOOST_TEST( string_to_unsigned( "12", v ) ) && BOOST_TEST_EQ( v, 12 );
    BOOST_TEST( string_to_unsigned( "123", v ) ) && BOOST_TEST_EQ( v, 123 );
    BOOST_TEST( string_to_unsigned( "1234", v ) ) && BOOST_TEST_EQ( v, 1234 );
    BOOST_TEST( string_to_unsigned( "12345", v ) ) && BOOST_TEST_EQ( v, 12345 );
    BOOST_TEST( string_to_unsigned( "123456", v ) ) && BOOST_TEST_EQ( v, 123456 );
    BOOST_TEST( string_to_unsigned( "1234567", v ) ) && BOOST_TEST_EQ( v, 1234567 );
    BOOST_TEST( string_to_unsigned( "12345678", v ) ) && BOOST_TEST_EQ( v, 12345678 );
    BOOST_TEST( string_to_unsigned( "123456789", v ) ) && BOOST_TEST_EQ( v, 123456789 );
    BOOST_TEST( string_to_unsigned( "1234567890", v ) ) && BOOST_TEST_EQ( v, 1234567890 );
    BOOST_TEST_NOT( string_to_unsigned( "12345678901", v ) ) && BOOST_TEST_EQ( v, 1234567890 );
    BOOST_TEST_NOT( string_to_unsigned( "123456789012", v ) ) && BOOST_TEST_EQ( v, 1234567890 );

    BOOST_TEST( string_to_unsigned( "4294967295", v ) ) && BOOST_TEST_EQ( v, 4294967295 );
    BOOST_TEST_NOT( string_to_unsigned( "4294967296", v ) ) && BOOST_TEST_EQ( v, 429496729 );

    BOOST_TEST( string_to_unsigned( "01", v ) ) && BOOST_TEST_EQ( v, 1 );
    BOOST_TEST( string_to_unsigned( "001", v ) ) && BOOST_TEST_EQ( v, 1 );
    BOOST_TEST( string_to_unsigned( "0001", v ) ) && BOOST_TEST_EQ( v, 1 );
    BOOST_TEST( string_to_unsigned( "00001", v ) ) && BOOST_TEST_EQ( v, 1 );
    BOOST_TEST( string_to_unsigned( "000001", v ) ) && BOOST_TEST_EQ( v, 1 );
    BOOST_TEST( string_to_unsigned( "0000001", v ) ) && BOOST_TEST_EQ( v, 1 );
    BOOST_TEST( string_to_unsigned( "00000001", v ) ) && BOOST_TEST_EQ( v, 1 );
    BOOST_TEST( string_to_unsigned( "000000001", v ) ) && BOOST_TEST_EQ( v, 1 );
    BOOST_TEST( string_to_unsigned( "0000000001", v ) ) && BOOST_TEST_EQ( v, 1 );
    BOOST_TEST( string_to_unsigned( "00000000001", v ) ) && BOOST_TEST_EQ( v, 1 );
    BOOST_TEST( string_to_unsigned( "000000000001", v ) ) && BOOST_TEST_EQ( v, 1 );
    BOOST_TEST( string_to_unsigned( "0000000000001", v ) ) && BOOST_TEST_EQ( v, 1 );
    BOOST_TEST( string_to_unsigned( "00000000000001", v ) ) && BOOST_TEST_EQ( v, 1 );

    BOOST_TEST_NOT( string_to_unsigned( "1a", v ) ) && BOOST_TEST_EQ( v, 1 );
    BOOST_TEST_NOT( string_to_unsigned( "2 ", v ) ) && BOOST_TEST_EQ( v, 2 );

    return boost::report_errors();
}
