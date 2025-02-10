// Boost.Geometry
// Unit Test

// Copyright (c) 2018-2019 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "test_buffer_geo.hpp"

static std::string const simplex = "POINT(4.9 52.0)";

template <typename Formula, bool Clockwise, typename PointType>
void test_point()
{
    using polygon = bg::model::polygon<PointType, Clockwise>;

    // Test with a high tolerance to account for possible differences in andoyer/thomas.
    auto make_settings = [](int points_per_circle)
    {
        ut_settings result;
        result.tolerance = 0.5;
        result.points_per_circle = points_per_circle;
        return result;
    };

    // Test with a high tolerance to account for possible differences in andoyer/thomas.
    auto make_circle = [](int points_per_circle)
    {
        return bg::strategy::buffer::geographic_point_circle<Formula>(points_per_circle);
    };

    bg::strategies::buffer::geographic<Formula> strategy;
    bg::strategy::buffer::geographic_join_miter<Formula> join;
    bg::strategy::buffer::geographic_side_straight<Formula> side;
    bg::strategy::buffer::end_flat end;

    test_one_geo<PointType, polygon>("simplex_1_16", simplex, strategy, side, make_circle(360), join, end, 3.1415, 1.0, make_settings(360));
    test_one_geo<PointType, polygon>("simplex_5_8", simplex, strategy, side, make_circle(8), join, end, 70.7107, 5.0, make_settings(8));
    test_one_geo<PointType, polygon>("simplex_5_16", simplex, strategy, side, make_circle(16), join, end, 76.5437, 5.0, make_settings(16));
    test_one_geo<PointType, polygon>("simplex_5_32", simplex, strategy, side, make_circle(32), join, end, 77.9640, 5.0, make_settings(32));

    // The more points used for the buffer, the more the area approaches 10*PI square meters
    test_one_geo<PointType, polygon>("simplex_10_8", simplex, strategy, side, make_circle(8), join, end, 282.8430, 10.0, make_settings(8));
    test_one_geo<PointType, polygon>("simplex_10_16", simplex, strategy, side, make_circle(16), join, end, 306.1471, 10.0, make_settings(16));
    test_one_geo<PointType, polygon>("simplex_10_32", simplex, strategy, side, make_circle(32), join, end, 312.1450, 10.0, make_settings(32));
    test_one_geo<PointType, polygon>("simplex_10_180", simplex, strategy, side, make_circle(180), join, end, 313.9051, 10.0, make_settings(180));
    test_one_geo<PointType, polygon>("simplex_10_180", simplex, strategy, side, make_circle(360), join, end, 314.15, 10.0, make_settings(360));
}

int test_main(int, char* [])
{
    BoostGeometryWriteTestConfiguration();

    test_point<bg::strategy::andoyer, true, bg::model::point<default_test_type, 2, bg::cs::geographic<bg::degree> > >();
    test_point<bg::strategy::thomas, true, bg::model::point<default_test_type, 2, bg::cs::geographic<bg::degree> > >();

#if ! defined(BOOST_GEOMETRY_TEST_ONLY_ONE_TYPE)
    test_point<bg::strategy::andoyer, true, bg::model::point<long double, 2, bg::cs::geographic<bg::degree> > >();
#endif

    return 0;
}
