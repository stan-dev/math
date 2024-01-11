// Copyright 2023 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/tuple/tuple.hpp>
#include <boost/core/lightweight_test.hpp>

int main()
{
    boost::tuple<int, int, int> tp( 1, 2, 3 );

    BOOST_TEST_EQ( boost::get<0>(tp), 1 );
    BOOST_TEST_EQ( boost::get<1>(tp), 2 );
    BOOST_TEST_EQ( boost::get<2>(tp), 3 );

    return boost::report_errors();
}
