// Copyright 2024 Ruben Perez Hidalgo.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/compat/to_array.hpp>
#include <boost/core/lightweight_test.hpp>

#include <array>
#include <vector>

namespace compat = boost::compat;

int main()
{
    {
        // regular C array
        int input[] = {1, 2};
        auto output = compat::to_array(input);
        static_assert(std::is_same<decltype(output), std::array<int, 2>>::value, "");
        BOOST_TEST_EQ(output[0], 1);
        BOOST_TEST_EQ(output[1], 2);
    }
    {
        // regular C array, const
        const int input[] = {1, 2};
        auto output = compat::to_array(input);
        static_assert(std::is_same<decltype(output), std::array<int, 2>>::value, "");
        BOOST_TEST_EQ(output[0], 1);
        BOOST_TEST_EQ(output[1], 2);
    }
    {
        // values not moved
        const std::vector<int> vec{1, 2};
        std::vector<int> input[]{vec};
        auto output = compat::to_array(input);
        static_assert(std::is_same<decltype(output), std::array<std::vector<int>, 1>>::value, "");
        BOOST_TEST_ALL_EQ(output[0].begin(), output[0].end(), vec.begin(), vec.end());
        BOOST_TEST_ALL_EQ(input[0].begin(), input[0].end(), vec.begin(), vec.end());  // input not modified
    }
    {
        // const values work
        const std::vector<int> vec{1, 2};
        const std::vector<int> input[]{vec};
        auto output = compat::to_array(input);
        static_assert(std::is_same<decltype(output), std::array<std::vector<int>, 1>>::value, "");
        BOOST_TEST_ALL_EQ(output[0].begin(), output[0].end(), vec.begin(), vec.end());
        BOOST_TEST_ALL_EQ(input[0].begin(), input[0].end(), vec.begin(), vec.end());  // input not modified
    }
    {
        // constexpr check
        constexpr int input[] = {1, 2, 3};
        constexpr auto output = compat::to_array(input);
        static_assert(std::is_same<decltype(output), const std::array<int, 3>>::value, "");
        BOOST_TEST_EQ(output[0], 1);
        BOOST_TEST_EQ(output[1], 2);
        BOOST_TEST_EQ(output[2], 3);
    }

    return boost::report_errors();
}
