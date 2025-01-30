// Boost.Geometry
// Unit Test Helper

// Copyright (c) 2018-2022 Barend Gehrels, Amsterdam, the Netherlands.

// This file was modified by Oracle on 2020-2022.
// Modifications copyright (c) 2020-2022 Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_GEOMETRY_TEST_BUFFER_GEO_HPP
#define BOOST_GEOMETRY_TEST_BUFFER_GEO_HPP

#include "test_buffer.hpp"

template<typename Geometry,
    typename GeometryOut,
    typename UmbrellaStrategy,
    typename SideStrategy,
    typename PointStrategy,
    typename JoinStrategy,
    typename EndStrategy>
void test_one_geo(std::string const& caseid,
        std::string const& wkt,
        UmbrellaStrategy const& strategy,
        SideStrategy const& side_strategy,
        PointStrategy const& circle_strategy,
        JoinStrategy const& join_strategy,
        EndStrategy const& end_strategy,
        double expected_area,
        double distance_left, ut_settings settings = ut_settings(),
        double distance_right = same_distance)
{
    Geometry input_geometry;
    bg::read_wkt(wkt, input_geometry);
    bg::correct(input_geometry);

    bool const symmetric = bg::math::equals(distance_right, same_distance);
    if (symmetric)
    {
        distance_right = distance_left;
    }
    auto const mean_distance = (distance_left + distance_right) / 2.0;

    bg::strategy::buffer::distance_asymmetric
    <
        typename bg::coordinate_type<Geometry>::type
    > distance_strategy(distance_left, distance_right);

    bg::model::multi_polygon<GeometryOut> buffer;

    test_buffer<GeometryOut>
            (caseid, buffer, input_geometry,
            join_strategy, end_strategy,
            distance_strategy, side_strategy, circle_strategy,
            strategy,
            -1, -1, expected_area,
            settings);

    if (symmetric && distance_left > 0.0)
    {
        // Verify if all the points of the output geometry are at or around the buffered distance
        // For linestrings with flat ends, it's not necessarily the case, there may be points
        // too close, especially on artefacts in heavily curved input with flat ends.
        // Therefore the default expectation can be modified. Inspect the SVG visually before doing this.
        std::size_t too_close = 0;
        std::size_t too_far = 0;
        std::size_t total = 0;
        boost::geometry::for_each_point(buffer, [&](const auto& p)
            {
                const auto distance = bg::distance(p, input_geometry);
                const auto f = distance / distance_left;

                if (f < 0.9) { too_close++; } else if (f > 1.1) { too_far++; }
                total++;
            });

        const double f = too_close / static_cast<double>(total);
        BOOST_CHECK_MESSAGE(f < settings.fraction_buffered_points_too_close,
                caseid << " has too many points too close " << too_close << " " << f);

        if (!JoinTestProperties<JoinStrategy>::is_miter())
        {
            BOOST_CHECK_MESSAGE(too_far == 0,
                    caseid << " has too far " << too_far);
        }
    }

    if (expected_area < 0 && bg::util::is_linear<Geometry>::value)
    {
        // Calculate the area of a linear feature using its length and the buffer distance.
        // For round ends, add the area of a circle (two halves at both ends).
        // For a straight line this expectation is perfect.
        // For a curved line it might be too large.
        // Therefore the default is 95% of it, and it can be modified with a setting.

        const auto addition = EndTestProperties<EndStrategy>::is_round()
            ? mean_distance * mean_distance * bg::math::pi<double>() : 0.0;
        const auto area = bg::area(buffer);
        const auto estimated_area = addition + bg::length(input_geometry) * (distance_left + distance_right);
        const auto min_area = settings.multiplier_min_area * estimated_area;
        const auto max_area = settings.multiplier_max_area * estimated_area;
        BOOST_CHECK_MESSAGE(area > min_area,
                caseid << " the area is too small, expected at least "
                << area << " " << min_area);
        BOOST_CHECK_MESSAGE(area < max_area,
                caseid << " the area is too large, expected at most "
                << area << " " << max_area);
    }
}

#endif
