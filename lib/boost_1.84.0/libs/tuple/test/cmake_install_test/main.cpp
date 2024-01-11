// Copyright 2017, 2021 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/tuple/tuple.hpp>
#include <cassert>

#define BOOST_TEST(expr) assert(expr)
#define BOOST_TEST_EQ(x1, x2) assert((x1)==(x2))

int main()
{
    boost::tuple<int, int, int> tp( 1, 2, 3 );

    BOOST_TEST_EQ( boost::get<0>(tp), 1 );
    BOOST_TEST_EQ( boost::get<1>(tp), 2 );
    BOOST_TEST_EQ( boost::get<2>(tp), 3 );
}
