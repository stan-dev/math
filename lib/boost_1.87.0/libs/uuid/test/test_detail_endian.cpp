// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/detail/endian.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cstdint>

int main()
{
    namespace detail = boost::uuids::detail;

    // byteswap u16

    {
        std::uint16_t x = 0x0102;
        std::uint16_t y = 0x0201;

        BOOST_TEST_EQ( detail::byteswap( x ), y );
        BOOST_TEST_EQ( detail::byteswap( y ), x );
    }

    {
        std::uint16_t x = 0xFFEE;
        std::uint16_t y = 0xEEFF;

        BOOST_TEST_EQ( detail::byteswap( x ), y );
        BOOST_TEST_EQ( detail::byteswap( y ), x );
    }

    // byteswap u32

    {
        std::uint32_t x = 0x01020304;
        std::uint32_t y = 0x04030201;

        BOOST_TEST_EQ( detail::byteswap( x ), y );
        BOOST_TEST_EQ( detail::byteswap( y ), x );
    }

    {
        std::uint32_t x = 0xFFEEDDCC;
        std::uint32_t y = 0xCCDDEEFF;

        BOOST_TEST_EQ( detail::byteswap( x ), y );
        BOOST_TEST_EQ( detail::byteswap( y ), x );
    }

    // byteswap u64

    {
        std::uint64_t x = 0x0102030405060708;
        std::uint64_t y = 0x0807060504030201;

        BOOST_TEST_EQ( detail::byteswap( x ), y );
        BOOST_TEST_EQ( detail::byteswap( y ), x );
    }

    {
        std::uint64_t x = 0xFFEEDDCCBBAA9988;
        std::uint64_t y = 0x8899AABBCCDDEEFF;

        BOOST_TEST_EQ( detail::byteswap( x ), y );
        BOOST_TEST_EQ( detail::byteswap( y ), x );
    }

    // byteswap u128

#if defined(__SIZEOF_INT128__)

    {
        __uint128_t x = ( static_cast<__uint128_t>( 0x0011223344556677 ) << 64 ) | 0x8899AABBCCDDEEFF;
        __uint128_t y = ( static_cast<__uint128_t>( 0xFFEEDDCCBBAA9988 ) << 64 ) | 0x7766554433221100;

        BOOST_TEST( detail::byteswap( x ) == y );
        BOOST_TEST( detail::byteswap( y ) == x );
    }

#endif

    // load u16

    {
        std::uint16_t x = 0x0102;
        std::uint16_t y = 0x0201;

        unsigned char data[] = { 0x01, 0x02 };

        BOOST_TEST_EQ( detail::load_little_u16( data ), y );
        BOOST_TEST_EQ( detail::load_big_u16( data ), x );

#if BOOST_UUID_BYTE_ORDER == BOOST_UUID_ORDER_LITTLE_ENDIAN

        BOOST_TEST_EQ( detail::load_native_u16( data ), y );

#else

        BOOST_TEST_EQ( detail::load_native_u16( data ), x );

#endif
    }

    {
        std::uint16_t x = 0xFFEE;
        std::uint16_t y = 0xEEFF;

        unsigned char data[] = { 0xFF, 0xEE };

        BOOST_TEST_EQ( detail::load_little_u16( data ), y );
        BOOST_TEST_EQ( detail::load_big_u16( data ), x );

#if BOOST_UUID_BYTE_ORDER == BOOST_UUID_ORDER_LITTLE_ENDIAN

        BOOST_TEST_EQ( detail::load_native_u16( data ), y );

#else

        BOOST_TEST_EQ( detail::load_native_u16( data ), x );

#endif
    }

    // load u32

    {
        std::uint32_t x = 0x01020304;
        std::uint32_t y = 0x04030201;

        unsigned char data[] = { 0x01, 0x02, 0x03, 0x04 };

        BOOST_TEST_EQ( detail::load_little_u32( data ), y );
        BOOST_TEST_EQ( detail::load_big_u32( data ), x );

#if BOOST_UUID_BYTE_ORDER == BOOST_UUID_ORDER_LITTLE_ENDIAN

        BOOST_TEST_EQ( detail::load_native_u32( data ), y );

#else

        BOOST_TEST_EQ( detail::load_native_u32( data ), x );

#endif
    }

    {
        std::uint32_t x = 0xFFEEDDCC;
        std::uint32_t y = 0xCCDDEEFF;

        unsigned char data[] = { 0xFF, 0xEE, 0xDD, 0xCC };

        BOOST_TEST_EQ( detail::load_little_u32( data ), y );
        BOOST_TEST_EQ( detail::load_big_u32( data ), x );

#if BOOST_UUID_BYTE_ORDER == BOOST_UUID_ORDER_LITTLE_ENDIAN

        BOOST_TEST_EQ( detail::load_native_u32( data ), y );

#else

        BOOST_TEST_EQ( detail::load_native_u32( data ), x );

#endif
    }

    // load u64

    {
        std::uint64_t x = 0x0102030405060708;
        std::uint64_t y = 0x0807060504030201;

        unsigned char data[] = { 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08 };

        BOOST_TEST_EQ( detail::load_little_u64( data ), y );
        BOOST_TEST_EQ( detail::load_big_u64( data ), x );

#if BOOST_UUID_BYTE_ORDER == BOOST_UUID_ORDER_LITTLE_ENDIAN

        BOOST_TEST_EQ( detail::load_native_u64( data ), y );

#else

        BOOST_TEST_EQ( detail::load_native_u64( data ), x );

#endif
    }

    {
        std::uint64_t x = 0xFFEEDDCCBBAA9988;
        std::uint64_t y = 0x8899AABBCCDDEEFF;

        unsigned char data[] = { 0xFF, 0xEE, 0xDD, 0xCC, 0xBB, 0xAA, 0x99, 0x88 };

        BOOST_TEST_EQ( detail::load_little_u64( data ), y );
        BOOST_TEST_EQ( detail::load_big_u64( data ), x );

#if BOOST_UUID_BYTE_ORDER == BOOST_UUID_ORDER_LITTLE_ENDIAN

        BOOST_TEST_EQ( detail::load_native_u64( data ), y );

#else

        BOOST_TEST_EQ( detail::load_native_u64( data ), x );

#endif
    }

    // load u128

#if defined(__SIZEOF_INT128__)

    {
        __uint128_t x = ( static_cast<__uint128_t>( 0x0011223344556677 ) << 64 ) | 0x8899AABBCCDDEEFF;
        __uint128_t y = ( static_cast<__uint128_t>( 0xFFEEDDCCBBAA9988 ) << 64 ) | 0x7766554433221100;

        unsigned char data[] = { 0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88, 0x99, 0xAA, 0xBB, 0xCC, 0xDD, 0xEE, 0xFF };

        BOOST_TEST( detail::load_little_u128( data ) == y );
        BOOST_TEST( detail::load_big_u128( data ) == x );

#if BOOST_UUID_BYTE_ORDER == BOOST_UUID_ORDER_LITTLE_ENDIAN

        BOOST_TEST( detail::load_native_u128( data ) == y );

#else

        BOOST_TEST( detail::load_native_u128( data ) == x );

#endif
    }

#endif

    // store u16

    {
        std::uint16_t x = 0x0102;
        std::uint16_t y = 0x0201;

        unsigned char data[ 2 ] = {};

        detail::store_little_u16( data, x );
        BOOST_TEST_EQ( detail::load_little_u16( data ), x );
        BOOST_TEST_EQ( detail::load_big_u16( data ), y );

        detail::store_big_u16( data, x );
        BOOST_TEST_EQ( detail::load_little_u16( data ), y );
        BOOST_TEST_EQ( detail::load_big_u16( data ), x );

        detail::store_native_u16( data, x );

#if BOOST_UUID_BYTE_ORDER == BOOST_UUID_ORDER_LITTLE_ENDIAN

        BOOST_TEST_EQ( detail::load_little_u16( data ), x );
        BOOST_TEST_EQ( detail::load_big_u16( data ), y );

#else

        BOOST_TEST_EQ( detail::load_little_u16( data ), y );
        BOOST_TEST_EQ( detail::load_big_u16( data ), x );

#endif
    }

    {
        std::uint16_t x = 0xFFEE;
        std::uint16_t y = 0xEEFF;

        unsigned char data[ 2 ] = {};

        detail::store_little_u16( data, x );
        BOOST_TEST_EQ( detail::load_little_u16( data ), x );
        BOOST_TEST_EQ( detail::load_big_u16( data ), y );

        detail::store_big_u16( data, x );
        BOOST_TEST_EQ( detail::load_little_u16( data ), y );
        BOOST_TEST_EQ( detail::load_big_u16( data ), x );

        detail::store_native_u16( data, x );

#if BOOST_UUID_BYTE_ORDER == BOOST_UUID_ORDER_LITTLE_ENDIAN

        BOOST_TEST_EQ( detail::load_little_u16( data ), x );
        BOOST_TEST_EQ( detail::load_big_u16( data ), y );

#else

        BOOST_TEST_EQ( detail::load_little_u16( data ), y );
        BOOST_TEST_EQ( detail::load_big_u16( data ), x );

#endif
    }

    // store u32

    {
        std::uint32_t x = 0x01020304;
        std::uint32_t y = 0x04030201;

        unsigned char data[ 4 ] = {};

        detail::store_little_u32( data, x );
        BOOST_TEST_EQ( detail::load_little_u32( data ), x );
        BOOST_TEST_EQ( detail::load_big_u32( data ), y );

        detail::store_big_u32( data, x );
        BOOST_TEST_EQ( detail::load_little_u32( data ), y );
        BOOST_TEST_EQ( detail::load_big_u32( data ), x );

        detail::store_native_u32( data, x );

#if BOOST_UUID_BYTE_ORDER == BOOST_UUID_ORDER_LITTLE_ENDIAN

        BOOST_TEST_EQ( detail::load_little_u32( data ), x );
        BOOST_TEST_EQ( detail::load_big_u32( data ), y );

#else

        BOOST_TEST_EQ( detail::load_little_u32( data ), y );
        BOOST_TEST_EQ( detail::load_big_u32( data ), x );

#endif
    }

    {
        std::uint32_t x = 0xFFEEDDCC;
        std::uint32_t y = 0xCCDDEEFF;

        unsigned char data[ 4 ] = {};

        detail::store_little_u32( data, x );
        BOOST_TEST_EQ( detail::load_little_u32( data ), x );
        BOOST_TEST_EQ( detail::load_big_u32( data ), y );

        detail::store_big_u32( data, x );
        BOOST_TEST_EQ( detail::load_little_u32( data ), y );
        BOOST_TEST_EQ( detail::load_big_u32( data ), x );

        detail::store_native_u32( data, x );

#if BOOST_UUID_BYTE_ORDER == BOOST_UUID_ORDER_LITTLE_ENDIAN

        BOOST_TEST_EQ( detail::load_little_u32( data ), x );
        BOOST_TEST_EQ( detail::load_big_u32( data ), y );

#else

        BOOST_TEST_EQ( detail::load_little_u32( data ), y );
        BOOST_TEST_EQ( detail::load_big_u32( data ), x );

#endif
    }

    // store u64

    {
        std::uint64_t x = 0x0102030405060708;
        std::uint64_t y = 0x0807060504030201;

        unsigned char data[ 8 ] = {};

        detail::store_little_u64( data, x );
        BOOST_TEST_EQ( detail::load_little_u64( data ), x );
        BOOST_TEST_EQ( detail::load_big_u64( data ), y );

        detail::store_big_u64( data, x );
        BOOST_TEST_EQ( detail::load_little_u64( data ), y );
        BOOST_TEST_EQ( detail::load_big_u64( data ), x );

        detail::store_native_u64( data, x );

#if BOOST_UUID_BYTE_ORDER == BOOST_UUID_ORDER_LITTLE_ENDIAN

        BOOST_TEST_EQ( detail::load_little_u64( data ), x );
        BOOST_TEST_EQ( detail::load_big_u64( data ), y );

#else

        BOOST_TEST_EQ( detail::load_little_u64( data ), y );
        BOOST_TEST_EQ( detail::load_big_u64( data ), x );

#endif
    }

    {
        std::uint64_t x = 0xFFEEDDCCBBAA9988;
        std::uint64_t y = 0x8899AABBCCDDEEFF;

        unsigned char data[ 8 ] = {};

        detail::store_little_u64( data, x );
        BOOST_TEST_EQ( detail::load_little_u64( data ), x );
        BOOST_TEST_EQ( detail::load_big_u64( data ), y );

        detail::store_big_u64( data, x );
        BOOST_TEST_EQ( detail::load_little_u64( data ), y );
        BOOST_TEST_EQ( detail::load_big_u64( data ), x );

        detail::store_native_u64( data, x );

#if BOOST_UUID_BYTE_ORDER == BOOST_UUID_ORDER_LITTLE_ENDIAN

        BOOST_TEST_EQ( detail::load_little_u64( data ), x );
        BOOST_TEST_EQ( detail::load_big_u64( data ), y );

#else

        BOOST_TEST_EQ( detail::load_little_u64( data ), y );
        BOOST_TEST_EQ( detail::load_big_u64( data ), x );

#endif
    }

    // store u128

#if defined(__SIZEOF_INT128__)

    {
        __uint128_t x = ( static_cast<__uint128_t>( 0x0011223344556677 ) << 64 ) | 0x8899AABBCCDDEEFF;
        __uint128_t y = ( static_cast<__uint128_t>( 0xFFEEDDCCBBAA9988 ) << 64 ) | 0x7766554433221100;

        unsigned char data[ 16 ] = {};

        detail::store_little_u128( data, x );
        BOOST_TEST( detail::load_little_u128( data ) == x );
        BOOST_TEST( detail::load_big_u128( data ) == y );

        detail::store_big_u128( data, x );
        BOOST_TEST( detail::load_little_u128( data ) == y );
        BOOST_TEST( detail::load_big_u128( data ) == x );

        detail::store_native_u128( data, x );

#if BOOST_UUID_BYTE_ORDER == BOOST_UUID_ORDER_LITTLE_ENDIAN

        BOOST_TEST( detail::load_little_u128( data ) == x );
        BOOST_TEST( detail::load_big_u128( data ) == y );

#else

        BOOST_TEST( detail::load_little_u128( data ) == y );
        BOOST_TEST( detail::load_big_u128( data ) == x );

#endif
    }

#endif

    return boost::report_errors();
}
