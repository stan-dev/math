// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/charconv.hpp>
#include <boost/core/lightweight_test.hpp>
#include <system_error>
#include <cstdint>
#include <cstring>

template <typename T>
void test()
{
    const char* str = "ffffffff1";

    for (int i = 11; i < 36; ++i)
    {
        T segment_result {42U};
        const auto r = boost::charconv::from_chars(str, str + std::strlen(str), segment_result, 16);
        BOOST_TEST(r.ec == std::errc::result_out_of_range);
        BOOST_TEST_EQ(segment_result, UINT32_C(42));
    }
}

template <typename T>
void test_64()
{
    const char* str = "ffffffffffffffff1";

    for (int i = 11; i < 36; ++i)
    {
        T segment_result {42U};
        const auto r = boost::charconv::from_chars(str, str + std::strlen(str), segment_result, 16);
        BOOST_TEST(r.ec == std::errc::result_out_of_range);
        BOOST_TEST_EQ(segment_result, UINT32_C(42));
    }
}

int main()
{
    test<std::uint32_t>();
    test<std::uint16_t>();
    test<std::int32_t>();
    test<std::int16_t>();

    test_64<std::uint64_t>();
    test_64<std::uint32_t>();
    test_64<std::int64_t>();
    test_64<std::int32_t>();

    return boost::report_errors();
}
