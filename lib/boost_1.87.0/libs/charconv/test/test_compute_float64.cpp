// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/charconv/detail/compute_float64.hpp>
#include <boost/core/lightweight_test.hpp>
#include <random>
#include <limits>
#include <cstdint>
#include <cmath>

using boost::charconv::detail::compute_float64;

inline void simple_test()
{
    bool success;

    // Trivial verifcation
    BOOST_TEST_EQ(compute_float64(1, 1, false, success), 1e1);
    BOOST_TEST_EQ(compute_float64(0, 1, true, success), -1e0);
    BOOST_TEST_EQ(compute_float64(308, 1, false, success), 1e308);

    // out of range
    BOOST_TEST_EQ(compute_float64(310, 5, false, success), HUGE_VAL);
    BOOST_TEST_EQ(compute_float64(-325, 5, false, success), 0);

    // Composite
    BOOST_TEST_EQ(compute_float64(10, 123456789, false, success), 123456789e10);
    BOOST_TEST_EQ(compute_float64(100, UINT64_C(4444444444444444444), false, success), 4444444444444444444e100);
}

int main(void)
{
    simple_test();

    return boost::report_errors();
}
