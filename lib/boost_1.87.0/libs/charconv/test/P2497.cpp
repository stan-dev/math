// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/charconv.hpp>
#include <boost/core/lightweight_test.hpp>
#include <system_error>
#include <type_traits>
#include <limits>
#include <cstring>
#include <cstdint>
#include <cerrno>
#include <utility>

// No overflows, negative numbers, locales, etc.
template <typename T>
void simple_from_chars_test()
{
    const char* buffer = "34";

    T v = 0;
    auto r = boost::charconv::from_chars(buffer, buffer + std::strlen(buffer), v);

    BOOST_TEST( r ) && BOOST_TEST_EQ(v, static_cast<T>(34));
    BOOST_TEST(r == r);

    boost::charconv::from_chars_result r2 {r.ptr, std::errc()};
    BOOST_TEST(r == r2);

    const char* buffer2 = "12";
    T v2 = 0;
    auto r3 = boost::charconv::from_chars(buffer2, buffer2 + std::strlen(buffer), v2);
    BOOST_TEST(r != r3);
    BOOST_TEST(r3) && BOOST_TEST_EQ(v2, static_cast<T>(12));
}

template <typename T>
void simple_to_chars_test()
{
    char buffer1[64] {};
    T v = 34;
    auto r1 = boost::charconv::to_chars(buffer1, buffer1 + sizeof(buffer1) - 1, v);
    BOOST_TEST(r1);
    BOOST_TEST_CSTR_EQ(buffer1, "34");

    boost::charconv::to_chars_result r {r1.ptr, r1.ec};
    BOOST_TEST(r1 == r);

    char buffer2[64] {};
    T v2 = 12;
    auto r2 = boost::charconv::to_chars(buffer2, buffer2 + sizeof(buffer2) - 1, v2);
    BOOST_TEST(r1 != r2);
    BOOST_TEST(r2);
    BOOST_TEST_CSTR_EQ(buffer2, "12");

    // If the base is not 2-36 inclusive return invalid value
    char buffer3[64] {};
    T v3 = 12;
    auto r3 = boost::charconv::to_chars(buffer3, buffer3 + sizeof(buffer3) - 1, v3, -2);
    BOOST_TEST(r3.ec == std::errc::invalid_argument);
    BOOST_TEST(!r3);
    auto r4 = boost::charconv::to_chars(buffer3, buffer3 + sizeof(buffer3) - 1, v3, 90);
    BOOST_TEST(r4.ec == std::errc::invalid_argument);
    BOOST_TEST(!r4);
}

int main()
{
    simple_from_chars_test<std::int32_t>();
    simple_from_chars_test<std::uint64_t>();

    simple_to_chars_test<std::int32_t>();
    simple_to_chars_test<std::uint64_t>();

    return boost::report_errors();
}

