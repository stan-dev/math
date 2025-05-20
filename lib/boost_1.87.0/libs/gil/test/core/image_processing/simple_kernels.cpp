//
// Copyright 2019 Olzhas Zhumabek <anonymous.from.applecity@gmail.com>
//
// Use, modification and distribution are subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
#include <boost/gil/image.hpp>
#include <boost/gil/image_processing/numeric.hpp>
#include <boost/gil/image_view.hpp>
#include <boost/gil/typedefs.hpp>

#include <boost/core/lightweight_test.hpp>

namespace gil = boost::gil;

void test_normalized_mean_generation()
{
    auto kernel = gil::generate_normalized_mean(5);
    for (const auto& cell: kernel)
    {
        const auto expected_value = static_cast<float>(1 / 25.f);
        BOOST_TEST_EQ(cell, expected_value);
    }
}

void test_unnormalized_mean_generation()
{
    auto kernel = gil::generate_unnormalized_mean(5);
    for (const auto& cell: kernel)
    {
        BOOST_TEST_EQ(cell, 1.0f);
    }
}

void test_gaussian_kernel_generation()
{
    auto kernel = boost::gil::generate_gaussian_kernel(7, 0.84089642);
    const float expected_values[7][7] =
    {
        {6.67847706e-07f, 2.29160778e-05f, 1.91169232e-04f, 3.87713181e-04f, 1.91169232e-04f, 2.29160778e-05f, 6.67847706e-07f},
        {2.29160778e-05f, 7.86326905e-04f, 6.55965268e-03f, 1.33037298e-02f, 6.55965268e-03f, 7.86326905e-04f, 2.29160778e-05f},
        {1.91169232e-04f, 6.55965268e-03f, 5.47215706e-02f, 1.10981636e-01f, 5.47215706e-02f, 6.55965268e-03f, 1.91169232e-04f},
        {3.87713181e-04f, 1.33037298e-02f, 1.10981636e-01f, 2.25083518e-01f, 1.10981636e-01f, 1.33037298e-02f, 3.87713181e-04f},
        {1.91169232e-04f, 6.55965268e-03f, 5.47215706e-02f, 1.10981636e-01f, 5.47215706e-02f, 6.55965268e-03f, 1.91169232e-04f},
        {2.29160778e-05f, 7.86326905e-04f, 6.55965268e-03f, 1.33037298e-02f, 6.55965268e-03f, 7.86326905e-04f, 2.29160778e-05f},
        {6.67847706e-07f, 2.29160778e-05f, 1.91169232e-04f, 3.87713181e-04f, 1.91169232e-04f, 2.29160778e-05f, 6.67847706e-07f}
    };

    double sum{0};
    for (gil::gray32f_view_t::coord_t y = 0; static_cast<std::size_t>(y) < kernel.size(); ++y)
    {
        for (gil::gray32f_view_t::coord_t x = 0; static_cast<std::size_t>(x) < kernel.size(); ++x)
        {
            auto output = kernel.at(static_cast<std::size_t>(x), static_cast<std::size_t>(y));
            auto expected = expected_values[y][x];
            auto percent_difference = std::ceil(std::abs(expected - output) / expected);
            BOOST_TEST_LT(percent_difference, 5);
            sum += output;
        }
    }
    BOOST_TEST_LT(std::abs(sum-1), 1e-7);
}

int main()
{
    test_normalized_mean_generation();
    test_unnormalized_mean_generation();
    test_gaussian_kernel_generation();
    return boost::report_errors();
}
