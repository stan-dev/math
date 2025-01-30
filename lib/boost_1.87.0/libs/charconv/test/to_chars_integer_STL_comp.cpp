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
#include <system_error>
#include <random>
#include <thread>
#include <vector>

template <typename T>
void stress_test_worker()
{
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<T> dist((std::numeric_limits<T>::min)(), (std::numeric_limits<T>::max)());

    std::mt19937_64 base_gen(rd());
    std::uniform_int_distribution<int> base_dist(2, 36);

    for (std::size_t i = 0; i < 100'000'000; ++i)
    {
        char buffer_stl[128] {};
        char buffer_boost[128] {};

        T v = dist(gen);
        int base = base_dist(base_gen);

        auto r_stl = std::to_chars(buffer_stl, buffer_stl + sizeof(buffer_stl) - 1, v, base);
        auto r_boost = boost::charconv::to_chars(buffer_boost, buffer_boost + sizeof(buffer_boost) - 1, v, base);

        BOOST_TEST_EQ(static_cast<std::ptrdiff_t>(r_stl.ptr - buffer_stl), static_cast<std::ptrdiff_t>(r_boost.ptr - buffer_boost));
        BOOST_TEST_CSTR_EQ(buffer_stl, buffer_boost);
    }
}

template <typename T>
void stress_test()
{
    auto num_threads = std::thread::hardware_concurrency();
    std::vector<std::thread> thread_manager(num_threads);

    for (std::size_t i = 0; i < num_threads; ++i)
    {
        thread_manager.emplace_back(std::thread(stress_test_worker<T>));
    }

    for (std::size_t i = 0; i < thread_manager.size(); ++i)
    {
        if (thread_manager[i].joinable())
        {
            thread_manager[i].join();
        }
    }
}

template <typename T, int base = 10>
void random_tests()
{
    std::mt19937_64 gen(42);
    std::uniform_int_distribution<T> dist((std::numeric_limits<T>::min)(), (std::numeric_limits<T>::max)());

    for (std::size_t i = 0; i < 100'000; ++i)
    {
        char buffer_stl[128] {};
        char buffer_boost[128] {};

        T v = dist(gen);

        auto r_stl = std::to_chars(buffer_stl, buffer_stl + sizeof(buffer_stl) - 1, v, base);
        auto r_boost = boost::charconv::to_chars(buffer_boost, buffer_boost + sizeof(buffer_boost) - 1, v, base);

        BOOST_TEST_EQ(static_cast<std::ptrdiff_t>(r_stl.ptr - buffer_stl), static_cast<std::ptrdiff_t>(r_boost.ptr - buffer_boost));
        BOOST_TEST_CSTR_EQ(buffer_stl, buffer_boost);
    }
}

template <typename T>
void test()
{
    // Signed 0
    char buffer1_stl[64] {};
    char buffer1_boost[64] {};
    T v1 = static_cast<T>(-0);
    auto r1_stl = std::to_chars(buffer1_stl, buffer1_stl + sizeof(buffer1_stl) - 1, v1);
    auto r1_boost = boost::charconv::to_chars(buffer1_boost, buffer1_boost + sizeof(buffer1_boost) - 1, v1);
    BOOST_TEST_EQ(r1_stl.ptr, buffer1_stl + 1);
    BOOST_TEST_EQ(r1_boost.ptr, buffer1_boost + 1);
    BOOST_TEST_CSTR_EQ(buffer1_stl, "0");
    BOOST_TEST_CSTR_EQ(buffer1_boost, "0");

    // Binary
    char buffer2_stl[64] {};
    char buffer2_boost[64] {};
    T v2 = static_cast<T>(42);
    auto r2_stl = std::to_chars(buffer2_stl, buffer2_stl + sizeof(buffer2_stl) - 1, v2, 2);
    auto r2_boost = boost::charconv::to_chars(buffer2_boost, buffer2_boost + sizeof(buffer2_boost) - 1, v2, 2);
    BOOST_TEST_EQ(r2_stl.ptr, buffer2_stl + 6);
    BOOST_TEST_EQ(r2_boost.ptr, buffer2_boost + 6);
    BOOST_TEST_CSTR_EQ(buffer2_stl, "101010");
    BOOST_TEST_CSTR_EQ(buffer2_boost, "101010");

    // Base 10
    char buffer3_stl[64] {};
    char buffer3_boost[64] {};
    T v3 = static_cast<T>(120);
    auto r3_stl = std::to_chars(buffer3_stl, buffer3_stl + sizeof(buffer3_stl) - 1, v3);
    auto r3_boost = boost::charconv::to_chars(buffer3_boost, buffer3_boost + sizeof(buffer3_boost) - 1, v3);
    BOOST_TEST_EQ(r3_stl.ptr, buffer3_stl + 3);
    BOOST_TEST_EQ(r3_boost.ptr, buffer3_boost + 3);
    BOOST_TEST_CSTR_EQ(buffer3_stl, "120");
    BOOST_TEST_CSTR_EQ(buffer3_boost, "120");

    // Hexadecimal
    char buffer4_stl[64] {};
    char buffer4_boost[64] {};
    T v4 = static_cast<T>(44);
    auto r4_stl = std::to_chars(buffer4_stl, buffer4_stl + sizeof(buffer4_stl) - 1, v4, 16);
    auto r4_boost = boost::charconv::to_chars(buffer4_boost, buffer4_boost + sizeof(buffer4_boost) - 1, v4, 16);
    BOOST_TEST_EQ(r4_stl.ptr, buffer4_stl + 2);
    BOOST_TEST_EQ(r4_boost.ptr, buffer4_boost + 2);
    BOOST_TEST_CSTR_EQ(buffer4_stl, "2c");
    BOOST_TEST_CSTR_EQ(buffer4_boost, "2c");

    // Base 28
    char buffer5_stl[64] {};
    char buffer5_boost[64] {};
    T v5 = static_cast<T>(100);
    auto r5_stl = std::to_chars(buffer5_stl, buffer5_stl + sizeof(buffer5_stl) - 1, v5, 28);
    auto r5_boost = boost::charconv::to_chars(buffer5_boost, buffer5_boost + sizeof(buffer5_boost) - 1, v5, 28);
    BOOST_TEST_EQ(r5_stl.ptr, buffer5_stl + 2);
    BOOST_TEST_EQ(r5_boost.ptr, buffer5_boost + 2);
    BOOST_TEST_CSTR_EQ(buffer5_stl, "3g");
    BOOST_TEST_CSTR_EQ(buffer5_boost, "3g");

    // Overflow
    char buffer6_stl[1] {};
    char buffer6_boost[1] {};
    T v6 = static_cast<T>(100);
    auto r6_stl = std::to_chars(buffer6_stl, buffer6_stl + sizeof(buffer6_stl) - 1, v6);
    auto r6_boost = boost::charconv::to_chars(buffer6_boost, buffer6_boost + sizeof(buffer6_boost) - 1, v6);
    BOOST_TEST_EQ(r6_stl.ptr, r6_stl.ptr);
    BOOST_TEST_EQ(r6_boost.ptr, r6_boost.ptr);
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

    // Specialized bases
    random_tests<int, 2>();
    random_tests<int, 4>();
    random_tests<int, 8>();
    random_tests<int, 16>();
    random_tests<int, 32>();

    #ifndef _MSC_VER // MSVC does not allow for 8-bit ints in std::uniform_int_distribution
    random_tests<std::int8_t, 10>();
    random_tests<std::uint8_t, 10>();
    #endif

    random_tests<std::int16_t, 10>();
    random_tests<std::uint16_t, 10>();
    random_tests<std::int32_t, 10>();
    random_tests<std::uint32_t, 10>();
    random_tests<std::int64_t, 10>();
    random_tests<std::uint64_t, 10>();

    // Generic implementation
    random_tests<int, 23>();

    // Stress tests are disabled for CI
    #ifdef BOOST_CHARCONV_STRESS_TEST
    stress_test<int>();
    #endif

    return boost::report_errors();
}

#else

int main()
{ 
    return 0;
}

#endif
