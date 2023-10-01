// Boost.Geometry
// Unit Test

// Copyright (c) 2022 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "test_buffer_geo.hpp"

// Testcase "simplex" is just simple (here non-closed) polygon.
static std::string const simplex = "POLYGON((10.4006587 63.4377982,10.4050904 63.4395993,10.4074994 63.4382527))";

// Testcase "Harbour" with several concave parts.
static std::string const harbour = "POLYGON((10.391199 63.436984,10.399588 63.440751,10.401838 63.441614,10.403849 63.444319,10.407099 63.443819,10.406201 63.441410,10.408428 63.441614,10.410235 63.443898,10.414224 63.444932,10.406201 63.435103,10.400019 63.435864,10.392565 63.435194,10.390452 63.434069,10.390258 63.435398,10.391156 63.435717,10.391474 63.436523,10.391838 63.435819,10.400974 63.439285,10.399258 63.440171,10.391319 63.436945))";

// Testcase "Nairobi" near the equator to verify with the QGis measurement tool if the buffer distance is right.
static std::string const nairobi = "POLYGON((36.8294584 -1.3193328,36.8298078 -1.3176255,36.8323390 -1.3181254,36.8320209 -1.3198413))";

template <typename Formula, bool Clockwise, typename PointType>
void test_geometry()
{
    using polygon = bg::model::polygon<PointType, Clockwise, true>;

    // Because areas can change significantly when another formula is used,
    // use a high tolerance.
    int const points_per_circle = 36;
    ut_settings settings(0.25);

    bg::strategies::buffer::geographic<Formula> strategy;
    bg::strategy::buffer::geographic_side_straight<Formula> side;
    bg::strategy::buffer::geographic_join_miter<Formula> join_miter;
    bg::strategy::buffer::geographic_join_round<Formula> join_round(points_per_circle);
    bg::strategy::buffer::geographic_point_circle<Formula> circle(points_per_circle);
    bg::strategy::buffer::end_flat end_flat;

    test_one_geo<polygon, polygon>("simplex_5_m", simplex, strategy, side, circle, join_miter, end_flat, 33015.0, 5.0, settings);
    test_one_geo<polygon, polygon>("simplex_5_r", simplex, strategy, side, circle, join_round, end_flat, 32940.0, 5.0, settings);

    test_one_geo<polygon, polygon>("harbour_5_m", harbour, strategy, side, circle, join_miter, end_flat, 545015.0, 5.0, settings);
    test_one_geo<polygon, polygon>("harbour_5_r", harbour, strategy, side, circle, join_round, end_flat, 544915.0, 5.0, settings);

    test_one_geo<polygon, polygon>("harbour_25_m", harbour, strategy, side, circle, join_miter, end_flat, 664343.0, 25.0, settings);
    test_one_geo<polygon, polygon>("harbour_25_r", harbour, strategy, side, circle, join_round, end_flat, 662205.0, 25.0, settings);

    test_one_geo<polygon, polygon>("harbour_50_m", harbour, strategy, side, circle, join_miter, end_flat, 818290.0, 50.0, settings);
    test_one_geo<polygon, polygon>("harbour_50_r", harbour, strategy, side, circle, join_round, end_flat, 805854.0, 50.0, settings);

    test_one_geo<polygon, polygon>("harbour_75_m", harbour, strategy, side, circle, join_miter, end_flat, 963064.0, 75.0, settings);
    test_one_geo<polygon, polygon>("harbour_75_r", harbour, strategy, side, circle, join_round, end_flat, 932937.0, 75.0, settings);

    test_one_geo<polygon, polygon>("harbour_30n_m", harbour, strategy, side, circle, join_miter, end_flat, 370066.0, -30.0, settings);
    test_one_geo<polygon, polygon>("harbour_30n_r", harbour, strategy, side, circle, join_round, end_flat, 371540.0, -30.0, settings);

    test_one_geo<polygon, polygon>("harbour_80n_m", harbour, strategy, side, circle, join_miter, end_flat, 185494.0, -80.0, settings);
    test_one_geo<polygon, polygon>("harbour_80n_r", harbour, strategy, side, circle, join_round, end_flat, 190110.0, -80.0, settings);

    test_one_geo<polygon, polygon>("nairobi_50_m", nairobi, strategy, side, circle, join_miter, end_flat, 113889.0, 50.0, settings);
    test_one_geo<polygon, polygon>("nairobi_50_r", nairobi, strategy, side, circle, join_round, end_flat, 111707.0, 50.0, settings);
}

int test_main(int, char* [])
{
    BoostGeometryWriteTestConfiguration();

    test_geometry<bg::strategy::andoyer, true, bg::model::point<default_test_type, 2, bg::cs::geographic<bg::degree> > >();

#if ! defined(BOOST_GEOMETRY_TEST_ONLY_ONE_TYPE)
    test_geometry<bg::strategy::thomas, true, bg::model::point<default_test_type, 2, bg::cs::geographic<bg::degree> > >();
    test_geometry<bg::strategy::andoyer, true, bg::model::point<long double, 2, bg::cs::geographic<bg::degree> > >();
#endif

    return 0;
}
