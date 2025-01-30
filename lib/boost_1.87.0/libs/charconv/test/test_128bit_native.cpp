// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/charconv/detail/emulated128.hpp>
#include <boost/core/lightweight_test.hpp>
#include <limits>
#include <cstdint>

void test128()
{
    auto r1 = boost::charconv::detail::umul128(1, 1);
    BOOST_TEST_EQ(r1.high, UINT64_C(0));
    BOOST_TEST_EQ(r1.low, UINT64_C(1));

    auto r2 = boost::charconv::detail::umul128(10, std::numeric_limits<std::uint64_t>::max());
    BOOST_TEST_EQ(r2.high, UINT64_C(9));
    BOOST_TEST_EQ(r2.low, UINT64_C(18446744073709551606));

    auto r3 = boost::charconv::detail::umul128(std::numeric_limits<std::uint64_t>::max(), std::numeric_limits<std::uint64_t>::max());
    BOOST_TEST_EQ(r3.high, UINT64_C(18446744073709551614));
    BOOST_TEST_EQ(r3.low, UINT64_C(1));
}

int main()
{
    test128();
    return boost::report_errors();
}
