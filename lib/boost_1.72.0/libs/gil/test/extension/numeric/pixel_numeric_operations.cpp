//
// Copyright 2019 Mateusz Loskot <mateusz at loskot dot net>
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//
#include <boost/gil.hpp>
#include <boost/gil/extension/numeric/pixel_numeric_operations.hpp>

#include <tuple>
#include <type_traits>

#define BOOST_TEST_MODULE test_ext_numeric_pixel_numeric_operations
#include "unit_test.hpp"
#include "unit_test_utility.hpp"
#include "core/image/test_fixture.hpp" // random_value
#include "core/pixel/test_fixture.hpp"

namespace gil = boost::gil;
namespace fixture = boost::gil::test::fixture;

BOOST_AUTO_TEST_SUITE(pixel_plus_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(plus_integer_same_types, pixel_t, fixture::pixel_integer_types)
{
    using channel_t = typename gil::channel_type<pixel_t>::type;
    gil::pixel_plus_t<pixel_t, pixel_t, pixel_t> f;
    {
        pixel_t p0;
        gil::static_fill(p0, static_cast<channel_t>(0));
        BOOST_TEST(f(p0, p0) == p0);
    }
    {
        pixel_t p1;
        gil::static_fill(p1, static_cast<channel_t>(1));
        pixel_t r2;
        gil::static_fill(r2, static_cast<channel_t>(2));
        BOOST_TEST(f(p1, p1) == r2);
    }
    {
        // Generates pixels with consecutive channel values: {1} or {1,2,3} or {1,2,3,4} etc.
        fixture::consecutive_value<channel_t> g(1);
        pixel_t p;
        gil::static_generate(p, [&g]() { return g(); });
        auto const r = f(p, p);
        BOOST_TEST(r != p);
        BOOST_TEST(gil::at_c<0>(r) == (gil::at_c<0>(p) + gil::at_c<0>(p)));
    }
}

BOOST_AUTO_TEST_SUITE_END() // pixel_plus_t

BOOST_AUTO_TEST_SUITE(pixel_minus_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(minus_integer_same_types, pixel_t, fixture::pixel_integer_types)
{
    using channel_t = typename gil::channel_type<pixel_t>::type;
    gil::pixel_minus_t<pixel_t, pixel_t, pixel_t> f;

    pixel_t p0;
    gil::static_fill(p0, static_cast<channel_t>(0));
    BOOST_TEST(f(p0, p0) == p0);
    {
        pixel_t p1, p2;
        gil::static_fill(p1, static_cast<channel_t>(1));
        gil::static_fill(p2, static_cast<channel_t>(2));
        pixel_t r1;
        gil::static_fill(r1, static_cast<channel_t>(1));
        BOOST_TEST(f(p2, p1) == r1);
    }
    {
        // Generates pixels with consecutive channel values: {1} or {1,2,3} or {1,2,3,4} etc.
        fixture::consecutive_value<channel_t> g(1);
        pixel_t p;
        gil::static_generate(p, [&g]() { return g(); });
        BOOST_TEST(f(p, p) == p0);
    }
}

BOOST_AUTO_TEST_SUITE_END() // pixel_minus_t

BOOST_AUTO_TEST_SUITE(pixel_multiplies_scalar_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(pixel_multiplies_scalar_integer_same_types, pixel_t, fixture::pixel_integer_types)
{
    using channel_t = typename gil::channel_type<pixel_t>::type;
    gil::pixel_multiplies_scalar_t<pixel_t, channel_t, pixel_t> f;

    pixel_t p0;
    gil::static_fill(p0, static_cast<channel_t>(0));
    BOOST_TEST(f(p0, 0) == p0);

    {
        pixel_t p1;
        gil::static_fill(p1, static_cast<channel_t>(1));
        BOOST_TEST(f(p1, 0) == p0);
        BOOST_TEST(f(p1, 1) == p1);
    }
    {
        // Generates pixels with consecutive channel values: {1} or {1,2,3} or {1,2,3,4} etc.
        fixture::consecutive_value<channel_t> g(1);
        pixel_t p;
        gil::static_generate(p, [&g]() { return g(); });

        // check first channel value is doubled
        auto const r = f(p, 2);
        BOOST_TEST(r != p);
        BOOST_TEST(gil::at_c<0>(r) == (gil::at_c<0>(p) * 2));
    }
}

BOOST_AUTO_TEST_SUITE_END() // pixel_multiplies_scalar_t

BOOST_AUTO_TEST_SUITE(pixel_multiply_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(pixel_multiply_integer_same_types, pixel_t, fixture::pixel_integer_types)
{
    using channel_t = typename gil::channel_type<pixel_t>::type;
    gil::pixel_multiply_t<pixel_t, pixel_t, pixel_t> f;

    pixel_t p0;
    gil::static_fill(p0, static_cast<channel_t>(0));
    BOOST_TEST(f(p0, p0) == p0);

    pixel_t p1;
    gil::static_fill(p1, static_cast<channel_t>(1));
    BOOST_TEST(f(p1, p1) == p1);

    pixel_t p2;
    gil::static_fill(p2, static_cast<channel_t>(2));
    BOOST_TEST(f(p1, p2) == p2);
}

BOOST_AUTO_TEST_SUITE_END() // pixel_multiply_t

BOOST_AUTO_TEST_SUITE(pixel_divides_scalar_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(pixel_divides_scalar_integer_same_types, pixel_t, fixture::pixel_integer_types)
{
    using channel_t = typename gil::channel_type<pixel_t>::type;
    gil::pixel_divides_scalar_t<pixel_t, channel_t, pixel_t> f;

    pixel_t p0;
    gil::static_fill(p0, static_cast<channel_t>(0));
    BOOST_TEST(f(p0, 1) == p0);

    pixel_t p1;
    gil::static_fill(p1, static_cast<channel_t>(1));
    BOOST_TEST(f(p1, 1) == p1);

    pixel_t p2;
    gil::static_fill(p2, static_cast<channel_t>(2));
    BOOST_TEST(f(p2, 2) == p1);
}

BOOST_AUTO_TEST_SUITE_END() // pixel_divides_scalar_t

BOOST_AUTO_TEST_SUITE(pixel_divide_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(pixel_divide_integer_same_types, pixel_t, fixture::pixel_integer_types)
{
    using channel_t = typename gil::channel_type<pixel_t>::type;
    gil::pixel_divide_t<pixel_t, pixel_t, pixel_t> f;

    pixel_t p0;
    gil::static_fill(p0, static_cast<channel_t>(0));
    pixel_t p1;
    gil::static_fill(p1, static_cast<channel_t>(1));
    BOOST_TEST(f(p0, p1) == p0);
    BOOST_TEST(f(p1, p1) == p1);

    pixel_t p2;
    gil::static_fill(p2, static_cast<channel_t>(2));
    BOOST_TEST(f(p2, p1) == p2);
}

BOOST_AUTO_TEST_SUITE_END() // pixel_divide_t

BOOST_AUTO_TEST_SUITE(pixel_halves_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(pixel_halves_integer_same_types, pixel_t, fixture::pixel_integer_types)
{
    using channel_t = typename gil::channel_type<pixel_t>::type;
    gil::pixel_halves_t<pixel_t> f;

    pixel_t p0;
    gil::static_fill(p0, static_cast<channel_t>(0));
    pixel_t p1;
    gil::static_fill(p1, static_cast<channel_t>(1));

    {
        auto p = p0;
        BOOST_TEST(f(p) == p0);
    }
    {
        auto p = p1;
        BOOST_TEST(f(p) == p0); // truncates toward Zero
    }
    {
        pixel_t p2;
        gil::static_fill(p2, static_cast<channel_t>(2));
        BOOST_TEST(f(p2) == p1);
    }
}

BOOST_AUTO_TEST_SUITE_END() // pixel_halves_t

BOOST_AUTO_TEST_SUITE(pixel_zeros_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(pixel_zeros_integer_same_types, pixel_t, fixture::pixel_integer_types)
{
    using channel_t = typename gil::channel_type<pixel_t>::type;
    gil::pixel_zeros_t<pixel_t> f;

    pixel_t p0;
    gil::static_fill(p0, static_cast<channel_t>(0));
    {
        auto p = p0;
        BOOST_TEST(f(p) == p0);
    }
    {
        fixture::consecutive_value<channel_t> g(1);
        pixel_t p;
        gil::static_generate(p, [&g]() { return g(); });
        BOOST_TEST(f(p) == p0);
    }
}

BOOST_AUTO_TEST_SUITE_END() // pixel_zeros_t

BOOST_AUTO_TEST_SUITE(zero_channels)

BOOST_AUTO_TEST_CASE_TEMPLATE(zero_channels_integer_same_types, pixel_t, fixture::pixel_integer_types)
{
    using channel_t = typename gil::channel_type<pixel_t>::type;

    pixel_t p0;
    gil::static_fill(p0, static_cast<channel_t>(0));
    {
        auto p = p0;
        gil::zero_channels(p);
        BOOST_TEST(p == p0);
    }
    {
        fixture::consecutive_value<channel_t> g(1);
        pixel_t p;
        gil::static_generate(p, [&g]() { return g(); });
        gil::zero_channels(p);
        BOOST_TEST(p == p0);
    }
}

BOOST_AUTO_TEST_SUITE_END() // zero_channels

BOOST_AUTO_TEST_SUITE(pixel_assigns_t)

BOOST_AUTO_TEST_CASE_TEMPLATE(pixel_assigns_integer_same_types, pixel_t, fixture::pixel_integer_types)
{
    using channel_t = typename gil::channel_type<pixel_t>::type;
    gil::pixel_assigns_t<pixel_t, pixel_t> f;

    {
        pixel_t p0, r;
        gil::static_fill(p0, static_cast<channel_t>(0));
        f(p0, r);
        BOOST_TEST(p0 == r);
    }
    {
        fixture::consecutive_value<channel_t> g(1);
        pixel_t p, r;
        gil::static_generate(p, [&g]() { return g(); });
        f(p, r);
        BOOST_TEST(p == r);
    }
}

BOOST_AUTO_TEST_SUITE_END() // pixel_assigns_t
