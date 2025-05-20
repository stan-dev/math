// Copyright 2023 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/system/result.hpp>
#include <boost/core/lightweight_test.hpp>

template<class T, class U> void test( T const& t, U const& u )
{
    BOOST_TEST_NE( static_cast<void const*>( &t ), static_cast<void const*>( &u ) );
}

int main()
{
    using namespace boost::system;

    test( result<int>::in_place_value, result<int>::in_place_error );
    test( result<void>::in_place_value, result<void>::in_place_error );
    test( result<int&>::in_place_value, result<int&>::in_place_error );

    return boost::report_errors();
}
