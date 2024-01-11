// Boost.Geometry (aka GGL, Generic Geometry Library)
// Unit Test

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.

// This file was modified by Oracle on 2021.
// Modifications copyright (c) 2021, Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_TEST_ENVELOPE_HPP
#define BOOST_GEOMETRY_TEST_ENVELOPE_HPP


#include <boost/variant/variant.hpp>

#include <geometry_test_common.hpp>

#include <boost/geometry/algorithms/envelope.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/geometry_collection.hpp>
#include <boost/geometry/strategies/strategies.hpp>

#include <boost/geometry/io/wkt/read.hpp>


template<typename Box, std::size_t DimensionCount>
struct check_result
{};

template <typename Box>
struct check_result<Box, 2>
{
    using ctype = typename bg::coordinate_type<Box>::type;
    using type = std::conditional_t
        <
            (std::is_integral<ctype>::value || std::is_floating_point<ctype>::value),
            double,
            ctype
        >;

    static void apply(Box const& b, type const& x1, type const& y1, type const& /*z1*/,
                type const& x2, type const& y2, type const& /*z2*/)
    {
        BOOST_CHECK_CLOSE((bg::get<bg::min_corner, 0>(b)), x1, 0.001);
        BOOST_CHECK_CLOSE((bg::get<bg::min_corner, 1>(b)), y1, 0.001);

        BOOST_CHECK_CLOSE((bg::get<bg::max_corner, 0>(b)), x2, 0.001);
        BOOST_CHECK_CLOSE((bg::get<bg::max_corner, 1>(b)), y2, 0.001);
    }
};

template <typename Box>
struct check_result<Box, 3>
{
    using ctype = typename bg::coordinate_type<Box>::type;
    using type = std::conditional_t
        <
            (std::is_integral<ctype>::value || std::is_floating_point<ctype>::value),
            double,
            ctype
        >;

    static void apply(Box const& b, type const& x1, type const& y1, type const& z1,
                type const& x2, type const& y2, type const& z2)
    {
        BOOST_CHECK_CLOSE((bg::get<bg::min_corner, 0>(b)), x1, 0.001);
        BOOST_CHECK_CLOSE((bg::get<bg::min_corner, 1>(b)), y1, 0.001);
        BOOST_CHECK_CLOSE((bg::get<bg::min_corner, 2>(b)), z1, 0.001);

        BOOST_CHECK_CLOSE((bg::get<bg::max_corner, 0>(b)), x2, 0.001);
        BOOST_CHECK_CLOSE((bg::get<bg::max_corner, 1>(b)), y2, 0.001);
        BOOST_CHECK_CLOSE((bg::get<bg::max_corner, 2>(b)), z2, 0.001);
    }
};


template <typename Geometry, typename T>
void test_envelope(std::string const& wkt,
                   T const& x1, T const& x2,
                   T const& y1, T const& y2,
                   T const& z1 = 0, T const& z2 = 0)
{
    typedef bg::model::box<typename bg::point_type<Geometry>::type > box_type;
    box_type b;

    Geometry geometry;
    bg::read_wkt(wkt, geometry);
    bg::envelope(geometry, b);
    check_result<box_type, bg::dimension<Geometry>::type::value>
            ::apply(b, x1, y1, z1, x2, y2, z2);

    boost::variant<Geometry> v(geometry);
    bg::envelope(v, b);
    check_result<box_type, bg::dimension<Geometry>::type::value>
            ::apply(b, x1, y1, z1, x2, y2, z2);

    bg::model::geometry_collection<boost::variant<Geometry>> gc{v};
    bg::envelope(gc, b);
    check_result<box_type, bg::dimension<Geometry>::type::value>
            ::apply(b, x1, y1, z1, x2, y2, z2);
}


#endif
