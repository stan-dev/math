//
// Copyright 2018 Mateusz Loskot <mateusz at loskot dot net>
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//
#include <boost/gil.hpp>

#include <boost/core/lightweight_test.hpp>

#include <cstdint>

namespace gil = boost::gil;

int main()
{
    gil::point_t d{2, 2};
    std::uint8_t r[] = { 1, 2, 3, 4 };
    std::uint8_t g[] = { 10, 20, 30, 40 };
    std::uint8_t b[] = { 110, 120, 130, 140 };
    std::uint8_t a[] = { 251, 252, 253, 254 };

    auto v = gil::planar_rgba_view(d.x, d.y, r, g, b, a, sizeof(std::uint8_t) * 2);
    BOOST_TEST(!v.empty());
    BOOST_TEST(v.dimensions() == d);
    BOOST_TEST(v.num_channels() == 4u);
    BOOST_TEST(v.size() == static_cast<std::size_t>(d.x * d.y));

    gil::rgba8_pixel_t const pf{1, 10, 110, 251};
    BOOST_TEST(v.front() == pf);

    gil::rgba8_pixel_t const pb{4, 40, 140, 254};
    BOOST_TEST(v.back() == pb);

    for (std::ptrdiff_t i = 0; i < v.size(); i++)
    {
        gil::rgba8_pixel_t const p{r[i], g[i], b[i], a[i]};
        BOOST_TEST(v[i] == p);
    }

    return boost::report_errors();
}
