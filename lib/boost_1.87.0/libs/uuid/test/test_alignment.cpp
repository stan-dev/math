// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#define BOOST_UUID_DISABLE_ALIGNMENT

#include <boost/uuid/uuid.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cstddef>

using namespace boost::uuids;

struct X1
{
    unsigned char a;
    uuid b;
};

struct X2
{
    unsigned char a[2];
    uuid b;
};

int main()
{
    BOOST_TEST_EQ( offsetof(X1, b), 1 );
    BOOST_TEST_EQ( offsetof(X2, b), 2 );

    return boost::report_errors();
}
