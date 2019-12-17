//
// Copyright 2005-2007 Adobe Systems Incorporated
// Copyright 2018 Mateusz Loskot <mateusz at loskot dot net>
//
// Distribtted under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//
#include <boost/gil/channel.hpp>
#include <boost/gil/channel_algorithm.hpp>
#include <boost/gil/typedefs.hpp>
#include <cstdint>
#include <limits>

#define BOOST_TEST_MODULE test_scoped_channel_value
#include "unit_test.hpp"

namespace gil = boost::gil;

struct int_minus_value  { static std::int8_t apply() { return -64; } };
struct int_plus_value   { static std::int8_t apply() { return  64; } };
using fixture = gil::scoped_channel_value
    <
        std::uint8_t, int_minus_value, int_plus_value
    >;

BOOST_AUTO_TEST_CASE(scoped_channel_value_default_constructor)
{
    fixture f;
    std::uint8_t v = f;
    BOOST_TEST(v == std::uint8_t{0});
}

BOOST_AUTO_TEST_CASE(scoped_channel_value_user_defined_constructors)
{
    fixture f{1};
    std::uint8_t v = f;
    BOOST_TEST(v == std::uint8_t{1});
}

BOOST_AUTO_TEST_CASE(scoped_channel_value_copy_constructors)
{
    fixture f1{128};
    fixture f2{f1};

    BOOST_TEST(std::uint8_t{f1} == std::uint8_t{128});
    BOOST_TEST(std::uint8_t{f1} == std::uint8_t{f2});
}

BOOST_AUTO_TEST_CASE(scoped_channel_value_assignment)
{
    fixture f;
    f = 64;
    std::uint8_t v = f;
    BOOST_TEST(v == std::uint8_t{64});
}

BOOST_AUTO_TEST_CASE(scoped_channel_value_float32_t)
{
    auto const tolerance = btt::tolerance(std::numeric_limits<float>::epsilon());
    // min
    BOOST_TEST(gil::float_point_zero<float>::apply() == 0.0, tolerance);
    BOOST_TEST(gil::channel_traits<gil::float32_t>::min_value() == 0.0);
    // max
    BOOST_TEST(gil::float_point_one<float>::apply() == 1.0, tolerance);
    BOOST_TEST(gil::channel_traits<gil::float32_t>::max_value() == 1.0);
}

BOOST_AUTO_TEST_CASE(scoped_channel_value_float64_t)
{
    auto const tolerance = btt::tolerance(std::numeric_limits<double>::epsilon());
    // min
    BOOST_TEST(gil::float_point_zero<double>::apply() == 0.0, tolerance);
    BOOST_TEST(gil::channel_traits<gil::float64_t>::min_value() == 0.0, tolerance);
    // max
    BOOST_TEST(gil::float_point_one<double>::apply() == 1.0, tolerance);
    BOOST_TEST(gil::channel_traits<gil::float64_t>::max_value() == 1.0, tolerance);
}

BOOST_AUTO_TEST_CASE(scoped_channel_value_halfs)
{
    // Create a double channel with range [-0.5 .. 0.5]
    struct minus_half { static double apply() { return -0.5; } };
    struct plus_half { static double apply() { return 0.5; } };
    using halfs = gil::scoped_channel_value<double, minus_half, plus_half>;

    auto const tolerance = btt::tolerance(std::numeric_limits<double>::epsilon());
    BOOST_TEST(gil::channel_traits<halfs>::min_value() == minus_half::apply(), tolerance);
    BOOST_TEST(gil::channel_traits<halfs>::max_value() == plus_half::apply(), tolerance);
    // scoped channel maximum should map to the maximum
    BOOST_TEST(gil::channel_convert<std::uint16_t>(
        gil::channel_traits<halfs>::max_value()) == 65535, tolerance);
}
