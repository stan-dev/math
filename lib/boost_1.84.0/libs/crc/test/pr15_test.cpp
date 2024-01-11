// Copyright 2023 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt)

#include <boost/crc.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cstring>

int main()
{
    char const* data = "123456789";

    {
        // CRC-14/DARC

        boost::crc_basic<14> c1( 0x805, 0x0000, 0x0000, true, true );
        boost::crc_optimal<14, 0x805, 0x0000, 0x0000, true, true> c2;

        c1.process_bytes( data, std::strlen( data ) );
        c2.process_bytes( data, std::strlen( data ) );

        BOOST_TEST_EQ( c1.checksum(), c2.checksum() );
    }

    {
        // CRC-24/BLE

        boost::crc_basic<24> c1( 0x00065B, 0x555555, 0x000000, true, true );
        boost::crc_optimal<24, 0x00065B, 0x555555, 0x000000, true, true> c2;

        c1.process_bytes( data, std::strlen( data ) );
        c2.process_bytes( data, std::strlen( data ) );

        BOOST_TEST_EQ( c1.checksum(), c2.checksum() );
    }

    return boost::report_errors();
}
