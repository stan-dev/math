// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/uuid.hpp>
#include <boost/config.hpp>

#if defined(BOOST_NO_CXX14_CONSTEXPR)

#include <boost/config/pragma_message.hpp>
BOOST_PRAGMA_MESSAGE("Skipping test because BOOST_NO_CXX14_CONSTEXPR is defined")

#else

#define STATIC_ASSERT(...) static_assert(__VA_ARGS__, #__VA_ARGS__)

using namespace boost::uuids;

constexpr uuid u1;

STATIC_ASSERT( u1.data()[ 0] == 0x00 );
STATIC_ASSERT( u1.data()[ 1] == 0x00 );
STATIC_ASSERT( u1.data()[ 2] == 0x00 );
STATIC_ASSERT( u1.data()[ 3] == 0x00 );
STATIC_ASSERT( u1.data()[ 4] == 0x00 );
STATIC_ASSERT( u1.data()[ 5] == 0x00 );
STATIC_ASSERT( u1.data()[ 6] == 0x00 );
STATIC_ASSERT( u1.data()[ 7] == 0x00 );
STATIC_ASSERT( u1.data()[ 8] == 0x00 );
STATIC_ASSERT( u1.data()[ 9] == 0x00 );
STATIC_ASSERT( u1.data()[10] == 0x00 );
STATIC_ASSERT( u1.data()[11] == 0x00 );
STATIC_ASSERT( u1.data()[12] == 0x00 );
STATIC_ASSERT( u1.data()[13] == 0x00 );
STATIC_ASSERT( u1.data()[14] == 0x00 );
STATIC_ASSERT( u1.data()[15] == 0x00 );

constexpr uuid u2 = {{ 0x01, 0x02, 0x03, 0x04 }};

STATIC_ASSERT( u2.data()[ 0] == 0x01 );
STATIC_ASSERT( u2.data()[ 1] == 0x02 );
STATIC_ASSERT( u2.data()[ 2] == 0x03 );
STATIC_ASSERT( u2.data()[ 3] == 0x04 );
STATIC_ASSERT( u2.data()[ 4] == 0x00 );
STATIC_ASSERT( u2.data()[ 5] == 0x00 );
STATIC_ASSERT( u2.data()[ 6] == 0x00 );
STATIC_ASSERT( u2.data()[ 7] == 0x00 );
STATIC_ASSERT( u2.data()[ 8] == 0x00 );
STATIC_ASSERT( u2.data()[ 9] == 0x00 );
STATIC_ASSERT( u2.data()[10] == 0x00 );
STATIC_ASSERT( u2.data()[11] == 0x00 );
STATIC_ASSERT( u2.data()[12] == 0x00 );
STATIC_ASSERT( u2.data()[13] == 0x00 );
STATIC_ASSERT( u2.data()[14] == 0x00 );
STATIC_ASSERT( u2.data()[15] == 0x00 );

constexpr uuid u3 = {{ 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F, 0x10 }};

STATIC_ASSERT( u3.data()[ 0] == 0x01 );
STATIC_ASSERT( u3.data()[ 1] == 0x02 );
STATIC_ASSERT( u3.data()[ 2] == 0x03 );
STATIC_ASSERT( u3.data()[ 3] == 0x04 );
STATIC_ASSERT( u3.data()[ 4] == 0x05 );
STATIC_ASSERT( u3.data()[ 5] == 0x06 );
STATIC_ASSERT( u3.data()[ 6] == 0x07 );
STATIC_ASSERT( u3.data()[ 7] == 0x08 );
STATIC_ASSERT( u3.data()[ 8] == 0x09 );
STATIC_ASSERT( u3.data()[ 9] == 0x0A );
STATIC_ASSERT( u3.data()[10] == 0x0B );
STATIC_ASSERT( u3.data()[11] == 0x0C );
STATIC_ASSERT( u3.data()[12] == 0x0D );
STATIC_ASSERT( u3.data()[13] == 0x0E );
STATIC_ASSERT( u3.data()[14] == 0x0F );
STATIC_ASSERT( u3.data()[15] == 0x10 );

#endif
