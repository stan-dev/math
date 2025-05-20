// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/charconv/detail/integer_search_trees.hpp>
#include <boost/charconv/detail/emulated128.hpp>
#include <boost/charconv/config.hpp>
#include <boost/core/lightweight_test.hpp>
#include <type_traits>
#include <limits>

#if defined(__GNUC__) && (__GNUC__ < 7)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Woverflow"
#endif

#ifdef BOOST_CHARCONV_HAS_INT128
// Test the u128 impl
template <typename T>
void test128()
{

    T max = BOOST_CHARCONV_UINT128_MAX;
    max /= 10;
    T val = 1;
    int i = 1;

    while (val < max)
    {
        BOOST_TEST_EQ(boost::charconv::detail::num_digits(val), i);
        ++i;
        val *= 10;
    }

    BOOST_TEST_EQ(boost::charconv::detail::num_digits(max * 10), 39);
}
#endif

#if defined(__GNUC__) && (__GNUC__ < 7)
# pragma GCC diagnostic pop
#endif

#if defined(__GNUC__) && (__GNUC__ >= 5)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wconversion"
#endif

// Tests the generic, u32, and u64 impls
template <typename T>
void test()
{
    constexpr T max = (std::numeric_limits<T>::max)() / 10;
    T val = 1;
    int i = 1;

    while (val < max)
    {
        BOOST_TEST_EQ(boost::charconv::detail::num_digits(val), i);
        ++i;
        val *= static_cast<T>(10);
    }

    BOOST_TEST_EQ(boost::charconv::detail::num_digits(val), std::numeric_limits<T>::digits10 + 1);
}

#if defined(__GNUC__) && (__GNUC__ >= 5)
# pragma GCC diagnostic pop
#endif

void test_emulated128()
{
    using namespace boost::charconv::detail;

    uint128 v1 {0, 1234};
    BOOST_TEST_EQ(num_digits(v1), 4);

    uint128 v2 {1234, 0};
    BOOST_TEST_EQ(num_digits(v2), 23);

    uint128 v3 {UINT64_MAX, UINT64_MAX};
    BOOST_TEST_EQ(num_digits(v3), 39);
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

    // Specialized tables
    test<std::uint32_t>();
    test<std::uint64_t>();
    #ifdef BOOST_CHARCONV_HAS_INT128
    test128<boost::uint128_type>();
    #endif

    test_emulated128();

    return boost::report_errors();
}
