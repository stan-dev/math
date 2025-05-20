// Boost.Geometry
// Unit Test

// Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <vector>

#include <boost/geometry/views/enumerate_view.hpp>
#include <boost/test/included/test_exec_monitor.hpp>

void test_const()
{
    const std::vector<std::string> vec{"cat", "mouse", "squirrel"};
    std::size_t test_index = 0;
    for (auto const& item : boost::geometry::util::enumerate(vec))
    {
        BOOST_CHECK_EQUAL(item.index, test_index++);
        switch(item.index)
        {
            case 0 : BOOST_CHECK_EQUAL(item.value, "cat"); break;
            case 1 : BOOST_CHECK_EQUAL(item.value, "mouse"); break;
            case 2 : BOOST_CHECK_EQUAL(item.value, "squirrel"); break;
        }
    }
}

void test_non_const()
{
    std::vector<std::string> vec{"Amsterdam", "London", "Paris"};
    std::size_t index_sum = 0;
    for (auto const& item : boost::geometry::util::enumerate(vec))
    {
        item.value += " is a city";
        index_sum += item.index;
    }
    BOOST_CHECK_EQUAL(vec[0], "Amsterdam is a city");
    BOOST_CHECK_EQUAL(vec[1], "London is a city");
    BOOST_CHECK_EQUAL(vec[2], "Paris is a city");
    BOOST_CHECK_EQUAL(index_sum, 3);
}

// Verifies the usage of the enumerate_view with C++17 structured bindings
// See https://en.cppreference.com/w/cpp/ranges/enumerate_view
void test_cpp17()
{
#if __cplusplus >= 201703L
    std::vector<int> numbers{1, 3, 5, 7};
    std::size_t sum_indexes = 0;
    int sum_numbers = 0;
    for (auto const [index, num] : boost::geometry::util::enumerate(numbers))
    {   
        sum_indexes += index;
        sum_numbers += num;
        // num is mutable even with const, which does not propagate to reference
        num++;
    }
    BOOST_CHECK_EQUAL(sum_indexes, 6);
    BOOST_CHECK_EQUAL(sum_numbers, 16);
    BOOST_CHECK_EQUAL(numbers[0], 2);
#endif
}

int test_main(int, char* [])
{
    test_const();
    test_non_const();
    test_cpp17();
    return 0;
}
