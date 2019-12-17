//
// Copyright 2019 Mateusz Loskot <mateusz at loskot dot net>
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//
#include <boost/gil/color_base.hpp>
#include <boost/gil/pixel.hpp>
#include <boost/gil/typedefs.hpp>

#include <type_traits>

#define BOOST_TEST_MODULE test_color_base_static_transform
#include "unit_test.hpp"

namespace gil = boost::gil;

BOOST_AUTO_TEST_CASE(single_source_gray8_to_gray8)
{
    gil::gray8_pixel_t src{128};
    gil::gray8_pixel_t dst{0};
    gil::static_transform(src, dst, [](std::uint8_t src_channel) {
        return src_channel; // copy
    });
    BOOST_TEST(gil::at_c<0>(src) == gil::at_c<0>(dst));
}

BOOST_AUTO_TEST_CASE(single_source_rgb8_to_rgb8)
{
    gil::rgb8_pixel_t src{32, 64, 128};
    gil::rgb8_pixel_t dst{0, 0, 0};
    gil::static_transform(src, dst, [](std::uint8_t src_channel) {
        return src_channel; // copy
    });
    BOOST_TEST(gil::at_c<0>(src) == gil::at_c<0>(dst));
    BOOST_TEST(gil::at_c<1>(src) == gil::at_c<1>(dst));
    BOOST_TEST(gil::at_c<2>(src) == gil::at_c<2>(dst));
}

BOOST_AUTO_TEST_CASE(single_source_rgb8_to_gray8)
{
    // Transformation of wider space to narrower space is a valid operation
    gil::rgb8_pixel_t  src{32,64, 128};
    gil::gray8_pixel_t dst{0};
    gil::static_transform(src, dst, [](std::uint8_t src_channel) {
        return src_channel; // copy
    });
    BOOST_TEST(gil::at_c<0>(dst) == std::uint8_t{32});
}

BOOST_AUTO_TEST_CASE(single_source_cmyk8_to_rgb8)
{
    // Transformation of wider space to narrower space is a valid operation
    gil::cmyk8_pixel_t src{16, 32, 64, 128};
    gil::rgb8_pixel_t dst{0, 0, 0};
    gil::static_transform(src, dst, [](std::uint8_t src_channel) {
        return src_channel; // copy
    });
    BOOST_TEST(gil::at_c<0>(dst) == std::uint8_t{16});
    BOOST_TEST(gil::at_c<1>(dst) == std::uint8_t{32});
    BOOST_TEST(gil::at_c<2>(dst) == std::uint8_t{64});
}

