// Copyright 2024 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/variant2/variant.hpp>
#include <boost/core/lightweight_test.hpp>

struct expr: boost::variant2::variant<bool, int, double>
{
    using variant::variant;
};

int main()
{
    using boost::variant2::get;

    {
        expr v{ true };

        BOOST_TEST_EQ( v.index(), 0 );
        BOOST_TEST_EQ( get<0>(v), true );
    }

    {
        expr v{ 2 };

        BOOST_TEST_EQ( v.index(), 1 );
        BOOST_TEST_EQ( get<1>(v), 2 );
    }

    {
        expr v{ 3.0 };

        BOOST_TEST_EQ( v.index(), 2 );
        BOOST_TEST_EQ( get<2>(v), 3.0 );
    }

    {
        expr v( true );

        BOOST_TEST_EQ( v.index(), 0 );
        BOOST_TEST_EQ( get<0>(v), true );
    }

    {
        expr v( 2 );

        BOOST_TEST_EQ( v.index(), 1 );
        BOOST_TEST_EQ( get<1>(v), 2 );
    }

    {
        expr v( 3.0 );

        BOOST_TEST_EQ( v.index(), 2 );
        BOOST_TEST_EQ( get<2>(v), 3.0 );
    }

    return boost::report_errors();
}
