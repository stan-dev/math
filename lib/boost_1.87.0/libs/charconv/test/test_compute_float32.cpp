// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/charconv/detail/compute_float32.hpp>
#include <boost/core/lightweight_test.hpp>
#include <limits>
#include <cstdint>
#include <cmath>

using boost::charconv::detail::compute_float32;

inline void simple_test()
{
    bool success;

    // Trivial verifcation
    BOOST_TEST_EQ(compute_float32(1, 1, false, success), 1e1F);
    BOOST_TEST_EQ(compute_float32(0, 1, true, success), -1e0F);
    BOOST_TEST_EQ(compute_float32(38, 1, false, success), 1e38F);

    // out of range
    BOOST_TEST_EQ(compute_float32(310, 5, false, success), HUGE_VALF);
    BOOST_TEST_EQ(compute_float32(-325, 5, false, success), 0.0F);

    // Composite
    BOOST_TEST_EQ(compute_float32(10, 123456789, false, success), 123456789e10F);
    BOOST_TEST_EQ(compute_float32(20, UINT64_C(444444444), false, success), 444444444e20F);
}

int main()
{
    simple_test();

    return boost::report_errors();
}
