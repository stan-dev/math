// Boost.Geometry
// Unit Test

// Copyright (c) 2022 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "test_buffer_geo.hpp"

// Case constructed such that an interior ring will be formed, depending on the buffer distance (11..18)
static std::string const trondheim = "MULTILINESTRING((10.3988246 63.4332944, 10.3991235 63.4341424, 10.3991956 63.4346470), (10.3985108 63.4336039, 10.3985342 63.4339071, 10.3984854 63.4342548, 10.3984154 63.4345156), (10.3973490 63.4332584, 10.4001666 63.4332774), (10.3985596 63.4346152, 10.3991129 63.4345792))";

static std::string const wrangel = "MULTILINESTRING((178.6083 71.0669,  179.9387 71.5239, -178.3429 71.5283, -177.4560 71.2289))";

template <typename Formula, bool Clockwise, typename PointType>
void test_geometry()
{
    using ml = bg::model::multi_linestring<bg::model::linestring<PointType>>;
    using polygon = bg::model::polygon<PointType, Clockwise>;

    int const points_per_circle = 360;
    ut_settings settings(0.5);
    settings.fraction_buffered_points_too_close = 0.35;

    bg::strategies::buffer::geographic<Formula> strategy;
    bg::strategy::buffer::geographic_side_straight<Formula> side;
    bg::strategy::buffer::geographic_join_miter<Formula> join_miter;
    bg::strategy::buffer::geographic_join_round<Formula> join_round(points_per_circle);
    bg::strategy::buffer::geographic_end_round<Formula> end_round(points_per_circle);
    bg::strategy::buffer::geographic_point_circle<Formula> circle(points_per_circle);
    bg::strategy::buffer::end_flat end_flat;

#if defined(BOOST_GEOMETRY_TEST_FAILURES)
    bool const andoyer_skip = false;
    bool const thomas_skip = false;
#else
    // TODO: some cases are missing one or more turns.
    bool const andoyer_skip = std::is_same<Formula, bg::strategy::andoyer>::value;
    bool const thomas_skip = std::is_same<Formula, bg::strategy::thomas>::value;
#endif

    test_one_geo<ml, polygon>("trondheim05_rr", trondheim, strategy, side, circle, join_round, end_round, 4398.0, 5.0, settings);

    if (! BOOST_GEOMETRY_CONDITION(thomas_skip))
    {
        test_one_geo<ml, polygon>("trondheim10_rr", trondheim, strategy, side, circle, join_round, end_round, 8994.0, 10.0, settings);
    }

#ifndef __APPLE__
    test_one_geo<ml, polygon>("trondheim12_rr", trondheim, strategy, side, circle, join_round, end_round, 10790.0, 12.0, settings);
#endif    

    if (! BOOST_GEOMETRY_CONDITION(thomas_skip) && ! BOOST_GEOMETRY_CONDITION(andoyer_skip))
    {
        test_one_geo<ml, polygon>("trondheim15_rr", trondheim, strategy, side, circle, join_round, end_round, 13358.0, 15.0, settings);
    }
    if (! BOOST_GEOMETRY_CONDITION(thomas_skip))
    {
        test_one_geo<ml, polygon>("trondheim17_rr", trondheim, strategy, side, circle, join_round, end_round, 14824.0, 17.0, settings);
    }

    test_one_geo<ml, polygon>("trondheim20_rr", trondheim, strategy, side, circle, join_round, end_round, 17055.0, 20.0, settings);
    test_one_geo<ml, polygon>("trondheim25_rr", trondheim, strategy, side, circle, join_round, end_round, 20657.0, 25.0, settings);

    test_one_geo<ml, polygon>("trondheim05_mf", trondheim, strategy, side, circle, join_miter, end_flat, 4190.0, 5.0, settings);
    test_one_geo<ml, polygon>("trondheim10_mf", trondheim, strategy, side, circle, join_miter, end_flat, 8196.0, 10.0, settings);
    test_one_geo<ml, polygon>("trondheim12_mf", trondheim, strategy, side, circle, join_miter, end_flat, 9706.0, 12.0, settings);
    test_one_geo<ml, polygon>("trondheim15_mf", trondheim, strategy, side, circle, join_miter, end_flat, 11708.0, 15.0, settings);
    test_one_geo<ml, polygon>("trondheim17_mf", trondheim, strategy, side, circle, join_miter, end_flat, 12896.0, 17.0, settings);
    test_one_geo<ml, polygon>("trondheim20_mf", trondheim, strategy, side, circle, join_miter, end_flat, 14486.0, 20.0, settings);
    test_one_geo<ml, polygon>("trondheim25_mf", trondheim, strategy, side, circle, join_miter, end_flat, 17025.0, 25.0, settings);
}

int test_main(int, char* [])
{
    BoostGeometryWriteTestConfiguration();

    test_geometry<bg::strategy::andoyer, true, bg::model::point<default_test_type, 2, bg::cs::geographic<bg::degree> > >();
    test_geometry<bg::strategy::thomas, true, bg::model::point<default_test_type, 2, bg::cs::geographic<bg::degree> > >();

#if ! defined(BOOST_GEOMETRY_TEST_ONLY_ONE_TYPE)
    test_geometry<true, bg::model::point<long double, 2, bg::cs::geographic<bg::degree> > >();
#endif

    return 0;
}
