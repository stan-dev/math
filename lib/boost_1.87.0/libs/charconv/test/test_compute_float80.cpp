// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/charconv/detail/compute_float80.hpp>
#include <boost/charconv/detail/emulated128.hpp>
#include <boost/charconv/detail/bit_layouts.hpp>
#include <boost/core/lightweight_test.hpp>
#include <random>
#include <limits>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <iomanip>

// MSVC uses long double = double
// Darwin sometimes uses double-double instead of long double
#if BOOST_CHARCONV_LDBL_BITS > 64 && !defined(__APPLE__) && !defined(_WIN32) && !defined(_WIN64)

using boost::charconv::detail::compute_float80;
using boost::charconv::detail::uint128;

template <typename T>
inline void test_fast_path()
{
    std::errc success;

    BOOST_TEST_EQ(compute_float80<long double>(1, T(1), false, success), 10.0L);
    BOOST_TEST_EQ(compute_float80<long double>(1, T(1), true, success), -10.0L);
    BOOST_TEST_EQ(compute_float80<long double>(27, T(1) << 112, false, success), 5.1922968585348276285304963292200960000000000000000e60L);
    BOOST_TEST_EQ(compute_float80<long double>(27, T(1) << 112, true, success), -5.1922968585348276285304963292200960000000000000000e60L);
}

int main()
{
    test_fast_path<uint128>();

    #ifdef BOOST_CHARCONV_HAS_INT128
    test_fast_path<boost::uint128_type>();
    #endif

    return boost::report_errors();
}

#else
int main()
{
    return 0;
}
#endif
