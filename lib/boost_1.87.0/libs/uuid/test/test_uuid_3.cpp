// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::uuids;

int main()
{
    {
        uuid u;

        std::uint8_t* p = u.data();

        BOOST_TEST_EQ( p, u.begin() );
        BOOST_TEST_EQ( p + u.size(), u.end() );
    }

    {
        uuid const u;

        std::uint8_t const* p = u.data();

        BOOST_TEST_EQ( p, u.begin() );
        BOOST_TEST_EQ( p + u.size(), u.end() );
    }

    return boost::report_errors();
}
