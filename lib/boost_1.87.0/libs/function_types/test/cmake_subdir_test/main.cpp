// Copyright 2017, 2021 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/function_types/function_arity.hpp>
#include <cassert>

#define BOOST_TEST(expr) assert(expr)
#define BOOST_TEST_EQ(x1, x2) assert((x1)==(x2))

int main()
{
    BOOST_TEST_EQ( boost::function_types::function_arity<void(int, int)>::value, 2 );
}
