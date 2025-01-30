// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/config.hpp>

#if !defined(BOOST_NO_CXX17_HDR_CHARCONV) && (!defined(__clang_major__) || (defined(__clang_major__) && __clang_major__ > 7))

#include <boost/charconv.hpp>
#include <boost/core/lightweight_test.hpp>
#include <charconv>
#include <type_traits>
#include <limits>
#include <cstring>
#include <cstdint>
#include <cerrno>
#include <utility>

template <typename T>
void test()
{
    // Binary
    const char* base2 = "0101010";
    T base2_v_boost = 0;
    T base2_v_stl = 0;
    auto r1_boost = boost::charconv::from_chars(base2, base2 + std::strlen(base2), base2_v_boost, 2);
    auto r1_stl = std::from_chars(base2, base2 + std::strlen(base2), base2_v_stl, 2);
    BOOST_TEST_EQ(r1_boost.ptr, r1_stl.ptr);
    BOOST_TEST_EQ(base2_v_boost, base2_v_stl);

    // Hexadecimal
    const char* base16 = "0x2a";
    T base16_v_boost = 0;
    T base16_v_stl = 0;
    auto r2_boost = boost::charconv::from_chars(base16, base16 + std::strlen(base16), base16_v_boost, 16);
    auto r2_stl = std::from_chars(base16, base16 + std::strlen(base16), base16_v_stl, 16);
    BOOST_TEST_EQ(r2_boost.ptr, r2_stl.ptr);
    BOOST_TEST_EQ(base16_v_boost, base16_v_stl);

    base16 = "2a";
    r2_boost = boost::charconv::from_chars(base16, base16 + std::strlen(base16), base16_v_boost, 16);
    r2_stl = std::from_chars(base16, base16 + std::strlen(base16), base16_v_stl, 16);
    BOOST_TEST_EQ(r2_boost.ptr, r2_stl.ptr);
    BOOST_TEST_EQ(base16_v_boost, base16_v_stl);

    // Decimal
    const char* base10 = "42";
    T base10_v_boost = 0;
    T base10_v_stl = 0;
    auto r3_boost = boost::charconv::from_chars(base10, base10 + std::strlen(base10), base10_v_boost);
    auto r3_stl = std::from_chars(base10, base10 + std::strlen(base10), base10_v_stl);
    BOOST_TEST_EQ(r3_boost.ptr, r3_stl.ptr);
    BOOST_TEST_EQ(base10_v_boost, base10_v_stl);

    // Negative number
    const char* neg = "-100";
    T negative_v_boost = 0;
    T negative_v_stl = 0;
    auto r4_boost = boost::charconv::from_chars(neg, neg + std::strlen(neg), negative_v_boost);
    auto r4_stl = std::from_chars(neg, neg + std::strlen(neg), negative_v_stl);
    BOOST_TEST_EQ(r4_boost.ptr, r4_stl.ptr);
    BOOST_TEST_EQ(negative_v_boost, negative_v_stl);

    // Digit seperator
    const char* sep = "1'000";
    T sep_v_boost = 0;
    T sep_v_stl = 0;
    auto r5_boost = boost::charconv::from_chars(sep, sep + std::strlen(sep), sep_v_boost);
    auto r5_stl = std::from_chars(sep, sep + std::strlen(sep), sep_v_stl);
    BOOST_TEST_EQ(r5_boost.ptr, r5_stl.ptr);
    BOOST_TEST_EQ(sep_v_boost, sep_v_stl);
}

int main()
{
    test<char>();
    test<signed char>();
    test<unsigned char>();
    test<short>();
    test<unsigned short>();
    test<int>();
    test<unsigned>();
    test<long>();
    test<unsigned long>();
    test<long long>();
    test<unsigned long long>();

    return boost::report_errors();
}

#else

int main()
{ 
    return 0;
}

#endif
