// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#define BOOST_UUID_DISABLE_ALIGNMENT

#include <boost/uuid/uuid.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/config/pragma_message.hpp>
#include <cstddef>

#if !defined(__GNUC__)

BOOST_PRAGMA_MESSAGE( "Test skipped, __GNUC__ is not defined" )
int main() {}

#else

using namespace boost::uuids;

// Note: recent GCCs produce the warning
// "ignoring packed attribute because of unpacked non-POD field"
// for X1 and X2, but this is a false positive; the packed attribute
// is not ignored and the tests pass.
//
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=114986

struct __attribute__((packed)) X1
{
    unsigned char a;
    uuid b;
    unsigned char c;
};

struct __attribute__((packed)) X2
{
    uuid a;
    unsigned char b;
    unsigned c;
    unsigned char d;
};

int main()
{
    BOOST_TEST_EQ( offsetof(X1, b), 1 );
    BOOST_TEST_EQ( sizeof(X1), 18 );

    BOOST_TEST_EQ( offsetof(X2, c), 17 );
    BOOST_TEST_EQ( sizeof(X2), 22 );

    return boost::report_errors();
}

#endif
