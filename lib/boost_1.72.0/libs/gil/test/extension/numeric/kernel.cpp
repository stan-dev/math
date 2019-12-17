//
// Copyright 2019 Mateusz Loskot <mateusz at loskot dot net>
// Copyright 2019 Miral Shah <miralshah2211@gmail.com>
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//
#define BOOST_DISABLE_ASSERTS 1 // kernel_1d_adaptor assertions are too strict
#include <boost/gil/extension/numeric/kernel.hpp>
#include <vector>

#define BOOST_TEST_MODULE test_ext_numeric_kernel
#include "unit_test.hpp"

namespace gil = boost::gil;

BOOST_AUTO_TEST_CASE(kernel_1d_default_constructor)
{
    gil::kernel_1d<int> k;
    BOOST_TEST(k.center() == 0);
    BOOST_TEST(k.left_size() == 0);
    BOOST_TEST(k.right_size() == -1); // TODO: Why not 0?
    // std::vector interface
    BOOST_TEST(k.size() == 0);
}

BOOST_AUTO_TEST_CASE(kernel_2d_default_constructor)
{
    gil::detail::kernel_2d<int> k;
    BOOST_TEST(k.center_y() == 0);
    BOOST_TEST(k.center_x() == 0);

    //BOOST_TEST(k.left_size() == 0);
    //BOOST_TEST(k.right_size() == -1);
    BOOST_TEST(k.upper_size() == 0);
    BOOST_TEST(k.lower_size() == -1);
    // std::vector interface
    BOOST_TEST(k.size() == 0);
}

BOOST_AUTO_TEST_CASE(kernel_1d_parameterized_constructor)
{
    gil::kernel_1d<int> k(9, 4);
    BOOST_TEST(k.center() == 4);
    BOOST_TEST(k.left_size() == 4);
    BOOST_TEST(k.right_size() == 4);
    // std::vector interface
    BOOST_TEST(k.size() == 9);
}

BOOST_AUTO_TEST_CASE(kernel_2d_parameterized_constructor)
{
    gil::detail::kernel_2d<int> k(9, 4, 4);
    BOOST_TEST(k.center_y() == 4);
    BOOST_TEST(k.center_x() == 4);
    BOOST_TEST(k.left_size() == 4);
    BOOST_TEST(k.right_size() == 4);
    BOOST_TEST(k.upper_size() == 4);
    BOOST_TEST(k.lower_size() == 4);
    // std::vector interface
    BOOST_TEST(k.size() == 9);
}

BOOST_AUTO_TEST_CASE(kernel_1d_parameterized_constructor_with_iterator)
{
    std::vector<int> v(9);
    gil::kernel_1d<int> k(v.cbegin(), v.size(), 4);
    BOOST_TEST(k.center() == 4);
    BOOST_TEST(k.left_size() == 4);
    BOOST_TEST(k.right_size() == 4);
    // std::vector interface
    BOOST_TEST(k.size() == 9);
}

BOOST_AUTO_TEST_CASE(kernel_2d_parameterized_constructor_with_iterator)
{
    std::vector<int> v(81);
    gil::detail::kernel_2d<int> k(v.cbegin(), v.size(), 4, 4);
    BOOST_TEST(k.center_y() == 4);
    BOOST_TEST(k.center_x() == 4);
    BOOST_TEST(k.left_size() == 4);
    BOOST_TEST(k.right_size() == 4);
    BOOST_TEST(k.upper_size() == 4);
    BOOST_TEST(k.lower_size() == 4);
    // std::vector interface
    BOOST_TEST(k.size() == 9);
}

BOOST_AUTO_TEST_CASE(kernel_1d_copy_constructor)
{
    gil::kernel_1d<int> d(9, 4);
    gil::kernel_1d<int> k(d);
    BOOST_TEST(k.center() == 4);
    BOOST_TEST(k.center() == d.center());
    BOOST_TEST(k.left_size() == d.left_size());
    BOOST_TEST(k.right_size() == d.right_size());
    // std::vector interface
    BOOST_TEST(k.size() == d.size());
}

BOOST_AUTO_TEST_CASE(kernel_2d_copy_constructor)
{
    gil::detail::kernel_2d<int> d(9, 4, 4);
    gil::detail::kernel_2d<int> k(d);
    BOOST_TEST(k.center_y() == 4);
    BOOST_TEST(k.center_x() == 4);
    BOOST_TEST(k.center_y() == d.center_y());
    BOOST_TEST(k.center_x() == d.center_x());
    BOOST_TEST(k.left_size() == d.left_size());
    BOOST_TEST(k.right_size() == d.right_size());
    BOOST_TEST(k.lower_size() == d.lower_size());
    BOOST_TEST(k.upper_size() == d.upper_size());
    // std::vector interface
    BOOST_TEST(k.size() == d.size());
}

BOOST_AUTO_TEST_CASE(kernel_1d_assignment_operator)
{
    gil::kernel_1d<int> d(9, 4);
    gil::kernel_1d<int> k;
    k = d;
    BOOST_TEST(k.center() == 4);
    BOOST_TEST(k.center() == d.center());
    BOOST_TEST(k.left_size() == d.left_size());
    BOOST_TEST(k.right_size() == d.right_size());
    // std::vector interface
    BOOST_TEST(k.size() == d.size());
}

BOOST_AUTO_TEST_CASE(kernel_2d_assignment_operator)
{
    gil::detail::kernel_2d<int> d(9, 4, 4);
    gil::detail::kernel_2d<int> k;
    k = d;
    BOOST_TEST(k.center_y() == 4);
    BOOST_TEST(k.center_x() == 4);
    BOOST_TEST(k.center_y() == d.center_y());
    BOOST_TEST(k.center_x() == d.center_x());
    BOOST_TEST(k.left_size() == d.left_size());
    BOOST_TEST(k.right_size() == d.right_size());
    BOOST_TEST(k.lower_size() == d.lower_size());
    BOOST_TEST(k.upper_size() == d.upper_size());
    // std::vector interface
    BOOST_TEST(k.size() == d.size());
}

BOOST_AUTO_TEST_CASE(kernel_1d_reverse_Kernel)
{
    gil::kernel_1d<int> d(12, 4);
    BOOST_TEST(d.center() == 4);
    auto k = gil::reverse_kernel(d);
    BOOST_TEST(k.center() == d.right_size());
    // std::vector interface
    BOOST_TEST(k.size() == d.size());
}

BOOST_AUTO_TEST_CASE(kernel_1d_fixed_default_constructor)
{
    gil::kernel_1d_fixed<int, 9> k;
    BOOST_TEST(k.center() == 0);
    BOOST_TEST(k.left_size() == 0);
    BOOST_TEST(k.right_size() == 8); // TODO: Why not 0 or -1 if not set?
    // std::array interface
    BOOST_TEST(k.size() == 9);
}

BOOST_AUTO_TEST_CASE(kernel_2d_fixed_default_constructor)
{
    gil::detail::kernel_2d_fixed<int, 9> k;
    BOOST_TEST(k.center_x() == 0);
    BOOST_TEST(k.center_y() == 0);
    BOOST_TEST(k.left_size() == 0);
    BOOST_TEST(k.right_size() == 8); // TODO: Why not 0 or -1 if not set?
    BOOST_TEST(k.upper_size() == 0);
    BOOST_TEST(k.lower_size() == 8);
    // std::array interface
    BOOST_TEST(k.size() == 9);
}

BOOST_AUTO_TEST_CASE(kernel_1d_fixed_parameterized_constructor)
{
    gil::kernel_1d_fixed<int, 9> k(4);
    BOOST_TEST(k.center() == 4);
    BOOST_TEST(k.left_size() == 4);
    BOOST_TEST(k.right_size() == 4);
    // std::vector interface
    BOOST_TEST(k.size() == 9);
}

BOOST_AUTO_TEST_CASE(kernel_2d_fixed_parameterized_constructor)
{
    gil::detail::kernel_2d_fixed<int, 9> k(4, 4);
    BOOST_TEST(k.center_x() == 4);
    BOOST_TEST(k.center_y() == 4);
    BOOST_TEST(k.left_size() == 4);
    BOOST_TEST(k.right_size() == 4);
    BOOST_TEST(k.upper_size() == 4);
    BOOST_TEST(k.lower_size() == 4);
    // std::vector interface
    BOOST_TEST(k.size() == 9);
}

BOOST_AUTO_TEST_CASE(kernel_1d_fixed_parameterized_constructor_with_iterator)
{
    // FIXME: The constructor should throw if v.size() < k.size()
    std::vector<int> v(9);
    gil::kernel_1d_fixed<int, 9> k(v.cbegin(), 4);
    BOOST_TEST((gil::kernel_1d_fixed<int, 9>::static_size) == 9);
    BOOST_TEST(k.center() == 4);
    BOOST_TEST(k.left_size() == 4);
    BOOST_TEST(k.right_size() == 4);
    // std::vector interface
    BOOST_TEST(k.size() == 9);
}

BOOST_AUTO_TEST_CASE(kernel_2d_fixed_parameterized_constructor_with_iterator)
{
//    // FIXME: The constructor should throw if v.size() < k.size()
    std::array<int, 81> v;
    gil::detail::kernel_2d_fixed<int, 9> k(v.cbegin(), 4, 4);
    BOOST_TEST((gil::detail::kernel_2d_fixed<int, 9>::static_size) == 9);
    BOOST_TEST(k.center_y() == 4);
    BOOST_TEST(k.center_x() == 4);
    BOOST_TEST(k.left_size() == 4);
    BOOST_TEST(k.right_size() == 4);
    BOOST_TEST(k.upper_size() == 4);
    BOOST_TEST(k.lower_size() == 4);
    // std::vector interface
    BOOST_TEST(k.size() == 9);
}

BOOST_AUTO_TEST_CASE(kernel_1d_fixed_copy_constructor)
{
    gil::kernel_1d_fixed<int, 9> d(4);
    gil::kernel_1d_fixed<int, 9> k(d);
    BOOST_TEST((gil::kernel_1d_fixed<int, 9>::static_size) == 9);
    BOOST_TEST(k.center() == 4);
    BOOST_TEST(k.center() == d.center());
    BOOST_TEST(k.left_size() == d.left_size());
    BOOST_TEST(k.right_size() == d.right_size());
    // std::vector interface
    BOOST_TEST(k.size() == d.size());
}

BOOST_AUTO_TEST_CASE(kernel_2d_fixed_copy_constructor)
{
    gil::detail::kernel_2d_fixed<int, 9> d(4, 4);
    gil::detail::kernel_2d_fixed<int, 9> k(d);
    BOOST_TEST((gil::detail::kernel_2d_fixed<int, 9>::static_size) == 9);
    BOOST_TEST(k.center_x() == 4);
    BOOST_TEST(k.center_y() == 4);
    BOOST_TEST(k.center_x() == d.center_x());
    BOOST_TEST(k.center_y() == d.center_y());
    BOOST_TEST(k.left_size() == d.left_size());
    BOOST_TEST(k.right_size() == d.right_size());
    BOOST_TEST(k.lower_size() == d.lower_size());
    BOOST_TEST(k.upper_size() == d.upper_size());
    // std::vector interface
    BOOST_TEST(k.size() == d.size());
}

BOOST_AUTO_TEST_CASE(kernel_1d_fixed_assignment_operator)
{
    gil::kernel_1d_fixed<int, 9> d(4);
    gil::kernel_1d_fixed<int, 9> k;
    k = d;
    BOOST_TEST((gil::kernel_1d_fixed<int, 9>::static_size) == 9);
    BOOST_TEST(k.center() == 4);
    BOOST_TEST(k.center() == d.center());
    BOOST_TEST(k.left_size() == d.left_size());
    BOOST_TEST(k.right_size() == d.right_size());
    // std::vector interface
    BOOST_TEST(k.size() == d.size());
}

BOOST_AUTO_TEST_CASE(kernel_2d_fixed_assignment_operator)
{
    gil::detail::kernel_2d_fixed<int, 9> d(4, 4);
    gil::detail::kernel_2d_fixed<int, 9> k;
    k = d;
    BOOST_TEST((gil::detail::kernel_2d_fixed<int, 9>::static_size) == 9);
    BOOST_TEST(k.center_x() == 4);
    BOOST_TEST(k.center_y() == 4);
    BOOST_TEST(k.center_x() == d.center_x());
    BOOST_TEST(k.center_y() == d.center_y());
    BOOST_TEST(k.left_size() == d.left_size());
    BOOST_TEST(k.right_size() == d.right_size());
    BOOST_TEST(k.lower_size() == d.lower_size());
    BOOST_TEST(k.upper_size() == d.upper_size());
    // std::vector interface
    BOOST_TEST(k.size() == d.size());
}

BOOST_AUTO_TEST_CASE(kernel_1d_fixed_reverse_Kernel)
{
    std::array<int, 3> values = {1, 2, 3};
    gil::kernel_1d_fixed<int, 3> d(values.begin(), 1);
    BOOST_TEST((gil::kernel_1d_fixed<int, 3>::static_size) == 3);
    BOOST_TEST(d == values);

    std::array<int, 3> values_rev = {3, 2, 1};
    auto k = gil::reverse_kernel(d);
    BOOST_TEST(k == values_rev);
}


