// Copyright 2023 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/system/result.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::system;

// Tricky mixed construction cases
// https://github.com/boostorg/system/issues/104
// https://brevzin.github.io//c++/2023/01/18/optional-construction/

int main()
{
    {
        result<int> r( make_error_code( errc::invalid_argument ) );
        result<result<int>> r2( r );

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( r2.value(), r );
    }

    {
        result<int> r( 5 );
        result<result<int>> r2( r );

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( r2.value(), r );
    }

    {
        result<int> const r( make_error_code( errc::invalid_argument ) );
        result<result<int>> r2( r );

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( r2.value(), r );
    }

    {
        result<int> const r( 5 );
        result<result<int>> r2( r );

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( r2.value(), r );
    }

    {
        result<int> r( make_error_code( errc::invalid_argument ) );
        result<result<int>> r2( std::move( r ) );

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( r2.value(), r );
    }

    {
        result<int> r( 5 );
        result<result<int>> r2( std::move( r ) );

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( r2.value(), r );
    }

    return boost::report_errors();
}
