//
// Copyright 2019 Mateusz Loskot <mateusz at loskot dot net>
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//
#include <boost/gil.hpp>
#include <boost/gil/extension/numeric/convolve.hpp>

#include <tuple>
#include <type_traits>

#define BOOST_TEST_MODULE test_ext_numeric_colvolve_rows
#include "unit_test.hpp"
#include "unit_test_utility.hpp"
#include "test_fixture.hpp"
#include "core/image/test_fixture.hpp"

namespace gil = boost::gil;
namespace fixture = boost::gil::test::fixture;

BOOST_AUTO_TEST_SUITE(convolve_rows)

BOOST_AUTO_TEST_CASE_TEMPLATE(image_1x1_kernel_1x1_identity, Image, fixture::image_types)
{
    auto const img = fixture::create_image<Image>(1, 1, 7);
    auto img_out = fixture::create_image<Image>(1, 1, 0);

    using pixel_t = typename Image::value_type;
    using channel_t = typename gil::channel_type<pixel_t>::type;
    auto const kernel = fixture::create_kernel<channel_t>({1});
    gil::convolve_rows<pixel_t>(const_view(img), kernel, view(img_out));

    // 1x1 kernel reduces convolution to multiplication
    BOOST_TEST(gil::const_view(img).front() == gil::const_view(img_out).front());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(image_1x1_kernel_3x3_identity, Image, fixture::image_types)
{
    auto const img = fixture::create_image<Image>(1, 1, 7);
    auto img_out = fixture::create_image<Image>(1, 1, 0);

    using pixel_t = typename Image::value_type;
    using channel_t = typename gil::channel_type<pixel_t>::type;
    auto const kernel = fixture::create_kernel<channel_t>({0, 0, 0, 0, 1, 0, 0, 0, 0});
    gil::convolve_rows<pixel_t>(const_view(img), kernel, view(img_out));

    BOOST_TEST(gil::const_view(img).front() == gil::const_view(img_out).front());
}

BOOST_AUTO_TEST_SUITE_END()
