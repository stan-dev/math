// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/uuid.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/config/pragma_message.hpp>
#include <cstddef>

#if !defined(__GNUC__) && !defined(_MSC_VER)

BOOST_PRAGMA_MESSAGE( "Test skipped, __GNUC__ and _MSC_VER are not defined" )
int main() {}

#else

using namespace boost::uuids;

#pragma pack( push, 1 )

struct X1
{
    unsigned char a;
    uuid b;
    unsigned char c;
};

struct X2
{
    uuid a;
    unsigned char b;
    unsigned c;
    unsigned char d;
};

#pragma pack( pop )

int main()
{
    BOOST_TEST_EQ( offsetof(X1, b), 1 );
    BOOST_TEST_EQ( sizeof(X1), 18 );

    BOOST_TEST_EQ( offsetof(X2, c), 17 );
    BOOST_TEST_EQ( sizeof(X2), 22 );

    return boost::report_errors();
}

#endif
