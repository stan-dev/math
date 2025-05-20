// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cstddef>
#include <cstring>
#include <cstdint>

using namespace boost::uuids;

int main()
{
    uuid u0 = {{}};

    BOOST_TEST( u0.is_nil() );
    BOOST_TEST_EQ( sizeof( u0.data ), 16u );

    uuid u1 = {{ 0x40, 0x41, 0x42, 0x43, 0x44, 0x45, 0x46, 0x47, 0x48, 0x49, 0x4A, 0x4B, 0x4C, 0x4D, 0x4E, 0x4F }};

    {
        uuid u2 = u0;
        BOOST_TEST_EQ( u2, u0 );

        std::memcpy( &u2, &u1, sizeof( u2 ) );
        BOOST_TEST_EQ( u2, u1 );
    }

    {
        uuid u2 = u0;
        BOOST_TEST_EQ( u2, u0 );

        std::memcpy( &u2.data, &u1.data, sizeof( u2.data ) );
        BOOST_TEST_EQ( u2, u1 );
    }

    {
        uuid u2 = u0;
        BOOST_TEST_EQ( u2, u0 );

        std::memcpy( u2.data, u1.data, sizeof( u2.data ) );
        BOOST_TEST_EQ( u2, u1 );
    }

    {
        uuid u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        for( int i = 0; i < 16; ++i )
        {
            std::uint8_t* p = u2.data + i;
            std::uint8_t& r = u2.data[ i ];

            BOOST_TEST_EQ( p, &r );
            BOOST_TEST_EQ( *p, 0x40 + i );
        }
    }

    {
        uuid u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        for( int i = 0; i < 16; ++i )
        {
            std::uint8_t const* p = u2.data + i;
            std::uint8_t const& r = u2.data[ i ];

            BOOST_TEST_EQ( p, &r );
            BOOST_TEST_EQ( *p, 0x40 + i );
        }
    }

    {
        uuid const u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        for( int i = 0; i < 16; ++i )
        {
            std::uint8_t const* p = u2.data + i;
            std::uint8_t const& r = u2.data[ i ];

            BOOST_TEST_EQ( p, &r );
            BOOST_TEST_EQ( *p, 0x40 + i );
        }
    }

    {
        uuid u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        unsigned char* p = reinterpret_cast<unsigned char*>( &u2.data );

        for( int i = 0; i < 16; ++i )
        {
            BOOST_TEST_EQ( p[ i ], 0x40 + i );
        }
    }

    {
        uuid u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        unsigned char const* p = reinterpret_cast<unsigned char const*>( &u2.data );

        for( int i = 0; i < 16; ++i )
        {
            BOOST_TEST_EQ( p[ i ], 0x40 + i );
        }
    }

    {
        uuid const u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        unsigned char const* p = reinterpret_cast<unsigned char const*>( &u2.data );

        for( int i = 0; i < 16; ++i )
        {
            BOOST_TEST_EQ( p[ i ], 0x40 + i );
        }
    }

    {
        uuid u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        unsigned char* p = reinterpret_cast<unsigned char*>( u2.data + 0 );

        for( int i = 0; i < 16; ++i )
        {
            BOOST_TEST_EQ( p[ i ], 0x40 + i );
        }
    }

    {
        uuid u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        unsigned char const* p = reinterpret_cast<unsigned char const*>( u2.data + 0 );

        for( int i = 0; i < 16; ++i )
        {
            BOOST_TEST_EQ( p[ i ], 0x40 + i );
        }
    }

    {
        uuid const u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        unsigned char const* p = reinterpret_cast<unsigned char const*>( u2.data + 0 );

        for( int i = 0; i < 16; ++i )
        {
            BOOST_TEST_EQ( p[ i ], 0x40 + i );
        }
    }

    {
        uuid u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        unsigned char* p = (unsigned char*)&u2.data;

        for( int i = 0; i < 16; ++i )
        {
            BOOST_TEST_EQ( p[ i ], 0x40 + i );
        }
    }

    {
        uuid u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        unsigned char const* p = (unsigned char const*)&u2.data;

        for( int i = 0; i < 16; ++i )
        {
            BOOST_TEST_EQ( p[ i ], 0x40 + i );
        }
    }

    {
        uuid const u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        unsigned char const* p = (unsigned char const*)&u2.data;

        for( int i = 0; i < 16; ++i )
        {
            BOOST_TEST_EQ( p[ i ], 0x40 + i );
        }
    }

    {
        uuid u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        unsigned char* p = (unsigned char*)u2.data;

        for( int i = 0; i < 16; ++i )
        {
            BOOST_TEST_EQ( p[ i ], 0x40 + i );
        }
    }

    {
        uuid u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        unsigned char const* p = (unsigned char const*)u2.data;

        for( int i = 0; i < 16; ++i )
        {
            BOOST_TEST_EQ( p[ i ], 0x40 + i );
        }
    }

    {
        uuid const u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        unsigned char const* p = (unsigned char const*)u2.data;

        for( int i = 0; i < 16; ++i )
        {
            BOOST_TEST_EQ( p[ i ], 0x40 + i );
        }
    }

    {
        uuid u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        void* p = &u2.data;

        BOOST_TEST_EQ( p, &u2 );
    }

    {
        uuid u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        void const* p = &u2.data;

        BOOST_TEST_EQ( p, &u2 );
    }

    {
        uuid const u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        void const* p = &u2.data;

        BOOST_TEST_EQ( p, &u2 );
    }

    {
        uuid u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        void* p = u2.data;

        BOOST_TEST_EQ( p, &u2 );
    }

    {
        uuid u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        void const* p = u2.data;

        BOOST_TEST_EQ( p, &u2 );
    }

    {
        uuid const u2 = u1;
        BOOST_TEST_EQ( u2, u1 );

        void const* p = u2.data;

        BOOST_TEST_EQ( p, &u2 );
    }

    return boost::report_errors();
}
