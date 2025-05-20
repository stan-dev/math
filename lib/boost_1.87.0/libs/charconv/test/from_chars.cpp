// Copyright 2022 Peter Dimov
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

#if defined(__has_include)
#  if __has_include(<string_view>)
#    include <string_view>
#  endif
#endif

#ifdef BOOST_CHARCONV_HAS_INT128
template <typename T>
void test_128bit_int()
{
    const char* buffer1 = "85070591730234615865843651857942052864"; // 2^126
    T test_value = 1;
    test_value = test_value << 126;
    T v1 = 0;
    auto r1 = boost::charconv::from_chars(buffer1, buffer1 + std::strlen(buffer1), v1);
    BOOST_TEST(r1.ec == std::errc());
    BOOST_TEST(v1 == test_value);

    #ifdef __GLIBCXX_TYPE_INT_N_0
    BOOST_TEST(std::numeric_limits<T>::max() > static_cast<T>(std::numeric_limits<unsigned long long>::max()));
    #else
    BOOST_IF_CONSTEXPR (std::is_same<T, boost::int128_type>::value)
    {
        BOOST_TEST(BOOST_CHARCONV_INT128_MAX > static_cast<T>(std::numeric_limits<unsigned long long>::max()));
    }
    else
    {
        BOOST_TEST(BOOST_CHARCONV_UINT128_MAX > static_cast<T>(std::numeric_limits<unsigned long long>::max()));
    }
    #endif
}

template <typename T>
void test_128bit_overflow();

template <>
void test_128bit_overflow<boost::int128_type>()
{
    const char* buffer1 = "170141183460469231731687303715884105728"; // max + 1
    boost::int128_type v1 = 1000;
    auto r1 = boost::charconv::from_chars(buffer1, buffer1 + std::strlen(buffer1), v1);
    BOOST_TEST(r1.ec == std::errc::result_out_of_range);

    const char* buffer2 = "-170141183460469231731687303715884105729"; // min - 1
    boost::int128_type v2 = 1000;
    auto r2 = boost::charconv::from_chars(buffer2, buffer2 + std::strlen(buffer2), v2);
    BOOST_TEST(r2.ec == std::errc::result_out_of_range);
}

template <>
void test_128bit_overflow<boost::uint128_type>()
{
    const char* buffer1 = "340282366920938463463374607431768211457"; // max + 1
    boost::uint128_type v1 = 1000;
    auto r1 = boost::charconv::from_chars(buffer1, buffer1 + std::strlen(buffer1), v1);
    BOOST_TEST(r1.ec == std::errc::result_out_of_range);

    const char* buffer2 = "-1"; // min - 1
    boost::uint128_type v2 = 1000;
    auto r2 = boost::charconv::from_chars(buffer2, buffer2 + std::strlen(buffer2), v2);
    BOOST_TEST(r2.ec == std::errc::invalid_argument);
}

#endif // 128-bit testing

#ifndef BOOST_NO_CXX14_CONSTEXPR
template <typename T>
constexpr std::pair<T, boost::charconv::from_chars_result> constexpr_test_helper()
{
    const char* buffer1 = "42";
    T v1 = 0;
    auto r1 = boost::charconv::from_chars(buffer1, buffer1 + 2, v1);

    return std::make_pair(v1, r1);
}

template <typename T>
constexpr void constexpr_test()
{
    constexpr auto results = constexpr_test_helper<T>();
    static_assert(results.second.ec == std::errc(), "No error");
    static_assert(results.first == 42, "Value is 42");
}

#endif

template <typename T>
void base2_test()
{
    // Includes leading 0 which should be ignored
    const char* buffer1 = "0101010";
    T v1 = 0;
    auto r1 = boost::charconv::from_chars(buffer1, buffer1 + std::strlen(buffer1), v1, 2);
    BOOST_TEST(r1.ec == std::errc());
    BOOST_TEST_EQ(v1, static_cast<T>(42));
}

template <typename T>
void base16_test()
{
    // In base 16 0x and 0X prefixes are not allowed
    const char* buffer1 = "2a";
    T v1 = 0;
    auto r1 = boost::charconv::from_chars(buffer1, buffer1 + std::strlen(buffer1), v1, 16);
    BOOST_TEST(r1.ec == std::errc());
    BOOST_TEST_EQ(v1, static_cast<T>(42));

    const char* buffer2 = "0";
    T v2 = 1;
    auto r2 = boost::charconv::from_chars(buffer2, buffer2 + std::strlen(buffer2), v2, 16);
    BOOST_TEST(r2.ec == std::errc());
    BOOST_TEST_EQ(v2, static_cast<T>(0));
}

template <typename T>
void overflow_test()
{
    const char* buffer1 = "1234";
    T v1 = 0;
    auto r1 = boost::charconv::from_chars(buffer1, buffer1 + std::strlen(buffer1), v1);

    BOOST_IF_CONSTEXPR((std::numeric_limits<T>::max)() < 1234)
    {
        BOOST_TEST(r1.ec == std::errc::result_out_of_range);
    }
    else
    {
        BOOST_TEST(r1.ec == std::errc()) && BOOST_TEST_EQ(v1, 1234);
    }

    const char* buffer2 = "123456789123456789123456789";
    T v2 = 0;
    auto r2 = boost::charconv::from_chars(buffer2, buffer2 + std::strlen(buffer2), v2);
    // In the event of overflow v2 is to be returned unmodified
    BOOST_TEST(r2.ec == std::errc::result_out_of_range) && BOOST_TEST_EQ(v2, static_cast<T>(0));
}

template <typename T>
void invalid_argument_test()
{
    const char* buffer1 = "";
    T v1 = 0;
    auto r1 = boost::charconv::from_chars(buffer1, buffer1 + std::strlen(buffer1), v1);
    BOOST_TEST(r1.ec == std::errc::invalid_argument);

    const char* buffer2 = "-";
    T v2 = 0;
    auto r2 = boost::charconv::from_chars(buffer2, buffer2 + std::strlen(buffer2), v2);
    BOOST_TEST(r2.ec == std::errc::invalid_argument);

    const char* buffer3 = "+";
    T v3 = 0;
    auto r3 = boost::charconv::from_chars(buffer3, buffer3 + std::strlen(buffer3), v3);
    BOOST_TEST(r3.ec == std::errc::invalid_argument);

    BOOST_IF_CONSTEXPR(std::is_unsigned<T>::value)
    {
        const char* buffer4 = "-123";
        T v4 = 0;
        auto r4 = boost::charconv::from_chars(buffer4, buffer4 + std::strlen(buffer4), v4);
        BOOST_TEST(r4.ec == std::errc::invalid_argument);
    }

    // Bases outside 2-36 inclusive return std::errc::invalid_argument
    const char* buffer5 = "23";
    T v5 = 0;
    auto r5 = boost::charconv::from_chars(buffer5, buffer5 + std::strlen(buffer5), v5, 1);
    BOOST_TEST(r5.ec == std::errc::invalid_argument);
    auto r6 = boost::charconv::from_chars(buffer5, buffer5 + std::strlen(buffer5), v5, 50);
    BOOST_TEST(r6.ec == std::errc::invalid_argument);

    const char* buffer7 = "+12345";
    T v7 = 3;
    auto r7 = boost::charconv::from_chars(buffer7, buffer7 + std::strlen(buffer7), v7);
    BOOST_TEST(r7.ec == std::errc::invalid_argument);
    BOOST_TEST_EQ(v7, static_cast<T>(3));

    const char* buffer8 = " 12345";
    T v8 = 3;
    auto r8 = boost::charconv::from_chars(buffer8, buffer8 + std::strlen(buffer8), v8);
    BOOST_TEST(r8.ec == std::errc::invalid_argument);
    BOOST_TEST_EQ(v8, static_cast<T>(3));

    const char* buffer9 = "123 45";
    T v9 = 3;
    auto r9 = boost::charconv::from_chars(buffer9, buffer9 + std::strlen(buffer9), v9);
    BOOST_TEST(r9);
    BOOST_TEST_EQ(v9, static_cast<T>(123));

    // https://github.com/boostorg/charconv/issues/188
    const char* buffer10 = "x";
    T o;
    auto r10 = boost::charconv::from_chars(buffer10, buffer10 + 1, o);
    BOOST_TEST(r10.ec == std::errc::invalid_argument);

    const char* buffer11 = "3x";
    auto r11 = boost::charconv::from_chars(buffer11, buffer11 + 2, o);
    BOOST_TEST(r11);
    BOOST_TEST_EQ(o, static_cast<T>(3));

    // https://github.com/boostorg/charconv/issues/213
    #if defined(__cpp_lib_string_view) && __cpp_lib_string_view >= 201606L
    T segment_result {12};
    std::string_view input{"/c"};
    auto r12 = boost::charconv::from_chars(input.data(), input.data() + input.size(), segment_result);
    BOOST_TEST(r12.ec == std::errc::invalid_argument);
    BOOST_TEST_EQ(segment_result, static_cast<T>(12));
    #endif
}

// No overflows, negative numbers, locales, etc.
template <typename T>
void simple_test()
{
    const char* buffer = "34";

    T v = 0;
    auto r = boost::charconv::from_chars(buffer, buffer + std::strlen(buffer), v);

    BOOST_TEST( r.ec == std::errc() ) && BOOST_TEST_EQ(v, static_cast<T>(34));
    BOOST_TEST(r == r);

    boost::charconv::from_chars_result r2 {r.ptr, std::errc()};
    BOOST_TEST(r == r2);

    const char* buffer2 = "12";
    T v2 = 0;
    auto r3 = boost::charconv::from_chars(buffer2, buffer2 + std::strlen(buffer), v2);
    BOOST_TEST(r != r3);
    BOOST_TEST(r3.ec == std::errc()) && BOOST_TEST_EQ(v2, static_cast<T>(12));
}

template <typename T>
void extended_ascii_codes()
{
    const char* buffer = "30±5"; // plus/minus is 177
    T v = 0;
    auto r = boost::charconv::from_chars(buffer, buffer + std::strlen(buffer), v);
    BOOST_TEST(r.ec == std::errc()) && BOOST_TEST_EQ(v, static_cast<T>(30));

    const char* buffer2 = "123°"; // Degrees is 186
    T v2 = 0;
    auto r2 = boost::charconv::from_chars(buffer2, buffer2 + std::strlen(buffer), v2);
    BOOST_TEST(r2.ec == std::errc()) && BOOST_TEST_EQ(v2, static_cast<T>(123));

    const char* buffer3 = "2¼"; // 1/4 is 188
    T v3 = 0;
    auto r3 = boost::charconv::from_chars(buffer3, buffer3 + std::strlen(buffer3), v3);
    BOOST_TEST(r3.ec == std::errc()) && BOOST_TEST_EQ(v3, static_cast<T>(2));

    const char* buffer4 = "123²"; // squared is 178
    T v4 = 0;
    auto r4 = boost::charconv::from_chars(buffer4, buffer4 + std::strlen(buffer4), v4);
    BOOST_TEST(r4.ec == std::errc()) && BOOST_TEST_EQ(v4, static_cast<T>(123));
}

int main()
{
    simple_test<char>();
    simple_test<signed char>();
    simple_test<unsigned char>();
    simple_test<short>();
    simple_test<unsigned short>();
    simple_test<int>();
    simple_test<unsigned>();
    simple_test<long>();
    simple_test<unsigned long>();
    simple_test<long long>();
    simple_test<unsigned long long>();
    simple_test<std::int32_t>();
    simple_test<std::uint64_t>();
    
    invalid_argument_test<int>();
    invalid_argument_test<unsigned>();
    invalid_argument_test<std::uint16_t>();

    overflow_test<char>();
    overflow_test<int>();

    base16_test<int>();
    base16_test<unsigned>();

    base2_test<unsigned char>();
    base2_test<long>();

    #if !(defined(__GNUC__) && __GNUC__ == 5)
    #   ifndef BOOST_NO_CXX14_CONSTEXPR
            constexpr_test<int>();
    #   endif
    #endif

    #ifdef BOOST_CHARCONV_HAS_INT128
    test_128bit_int<boost::int128_type>();
    test_128bit_int<boost::uint128_type>();

    test_128bit_overflow<boost::int128_type>();
    test_128bit_overflow<boost::uint128_type>();
    #endif

    extended_ascii_codes<int>();
    extended_ascii_codes<unsigned>();
    extended_ascii_codes<char>();
    extended_ascii_codes<unsigned char>();

    return boost::report_errors();
}
