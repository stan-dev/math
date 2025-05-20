// Copyright 2024 Ruben Perez Hidalgo.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/compat/to_array.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/config.hpp>
#include <boost/config/workaround.hpp>
#include <boost/config/pragma_message.hpp>

#include <array>
#include <memory>
#include <utility>
#include <vector>

#if defined(BOOST_MSVC) && BOOST_MSVC < 1910

BOOST_PRAGMA_MESSAGE( "Test skipped because BOOST_MSVC < 1910" )
int main() {}

#else

namespace compat = boost::compat;

int main()
{
    {
        // regular C array
        int input[] = {5, 6};
        auto output = compat::to_array(std::move(input));
        static_assert(std::is_same<decltype(output), std::array<int, 2>>::value, "");
        BOOST_TEST_EQ(output[0], 5);
        BOOST_TEST_EQ(output[1], 6);
    }
    {
        // regular C array, const
        const int input[] = {5, 6};
        auto output = compat::to_array(std::move(input));
        static_assert(std::is_same<decltype(output), std::array<int, 2>>::value, "");
        BOOST_TEST_EQ(output[0], 5);
        BOOST_TEST_EQ(output[1], 6);
    }
    {
        // values moved
        const std::vector<int> vec{1, 2};
        std::vector<int> input[]{vec};
        auto output = compat::to_array(std::move(input));
        static_assert(std::is_same<decltype(output), std::array<std::vector<int>, 1>>::value, "");
        BOOST_TEST_ALL_EQ(output[0].begin(), output[0].end(), vec.begin(), vec.end());
#if !BOOST_WORKAROUND(BOOST_MSVC, < 1920)
        BOOST_TEST(input[0].empty());  // input left in a moved-from state
#endif
    }
    {
        // const values work (although they don't result in an effective move)
        const std::vector<int> vec{1, 2};
        const std::vector<int> input[]{vec};
        auto output = compat::to_array(std::move(input));
        static_assert(std::is_same<decltype(output), std::array<std::vector<int>, 1>>::value, "");
        BOOST_TEST_ALL_EQ(output[0].begin(), output[0].end(), vec.begin(), vec.end());
        BOOST_TEST_ALL_EQ(input[0].begin(), input[0].end(), vec.begin(), vec.end());  // input not modified
    }
#if !BOOST_WORKAROUND(BOOST_MSVC, < 1920)
    {
        // move-only types work
        std::unique_ptr<int> input[] = {std::unique_ptr<int>{new int(42)}};
        int* ptr = input[0].get();
        auto output = compat::to_array(std::move(input));
        static_assert(std::is_same<decltype(output), std::array<std::unique_ptr<int>, 1>>::value, "");
        BOOST_TEST_EQ(output[0].get(), ptr);
        BOOST_TEST_EQ(input[0].get(), nullptr);  // input left in a moved-from state
    }
#endif

    {
        // constexpr check
        constexpr int input[] = {1, 2, 3};
        constexpr auto output = compat::to_array(std::move(input));
        static_assert(std::is_same<decltype(output), const std::array<int, 3>>::value, "");
        BOOST_TEST_EQ(output[0], 1);
        BOOST_TEST_EQ(output[1], 2);
        BOOST_TEST_EQ(output[2], 3);
    }

    return boost::report_errors();
}

#endif
