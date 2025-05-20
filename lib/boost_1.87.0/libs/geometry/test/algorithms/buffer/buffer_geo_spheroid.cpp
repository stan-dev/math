// Boost.Geometry
// Unit Test

// Copyright (c) 2022 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2023, Oracle and/or its affiliates.
// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "test_buffer_geo.hpp"

namespace
{
    int const points_per_circle = 360;
    std::string const road = "LINESTRING(10.3966569 63.4276957,10.3998059 63.4279182,10.4000859 63.4283889,10.3982915 63.4284015,10.3980902 63.4288772,10.3987772 63.4288520)";
    std::string const torg = "POINT(10.3937759 63.4302323)";
}

template
<
    typename FormulaPolicy,
    typename Spheroid,
    typename CalculationType
>
struct geo_buffer_accurate_area
    : public bg::strategies::buffer::geographic<FormulaPolicy, Spheroid, CalculationType>
{
    using base_t = bg::strategies::buffer::geographic<FormulaPolicy, Spheroid, CalculationType>;

public:

    geo_buffer_accurate_area() = default;

    explicit geo_buffer_accurate_area(Spheroid const& spheroid)
        : base_t(spheroid)
    {}

    template <typename Geometry>
    auto area(Geometry const&,
                std::enable_if_t<! bg::util::is_box<Geometry>::value> * = nullptr) const
    {
        return bg::strategy::area::geographic
            <
                bg::strategy::karney,
                bg::strategy::default_order<bg::strategy::karney>::value,
                Spheroid, CalculationType
            >(base_t::m_spheroid);
    }
};

template <typename Formula, bool Clockwise, typename Point, typename Spheroid>
void test_linestring(std::string const& label, Spheroid const& spheroid,
                     double expected_area_round, double expected_area_miter)
{
    using linestring = bg::model::linestring<Point>;
    using polygon = bg::model::polygon<Point, Clockwise>;

#ifdef __APPLE__
    ut_settings settings(1.0);
#else    
    ut_settings settings(0.1);
#endif    

    using CT = typename bg::coordinate_type<Point>::type;
    geo_buffer_accurate_area<Formula, Spheroid, CT> strategy(spheroid);

    bg::strategy::buffer::geographic_side_straight<Formula, Spheroid> side(spheroid);
    bg::strategy::buffer::geographic_join_miter<Formula, Spheroid> join_miter(spheroid);
    bg::strategy::buffer::geographic_join_round<Formula, Spheroid> join_round(spheroid, points_per_circle);
    bg::strategy::buffer::geographic_end_round<Formula, Spheroid> end_round(spheroid, points_per_circle);

    // Ignored for linear or areal features
    bg::strategy::buffer::geographic_point_circle<Formula, Spheroid> circle(spheroid, points_per_circle);

    test_one_geo<linestring, polygon>(label + "_round", road, strategy, side, circle, join_round, end_round, expected_area_round, 10.0, settings);
    test_one_geo<linestring, polygon>(label + "_miter", road, strategy, side, circle, join_miter, end_round, expected_area_miter, 10.0, settings);
}

template <typename Formula, bool Clockwise, typename Point, typename Spheroid>
void test_point(std::string const& label, Spheroid const& spheroid, double expected_area)
{
    using polygon = bg::model::polygon<Point, Clockwise>;

#ifdef __APPLE__
    ut_settings settings(1.15);
#else    
    ut_settings settings(0.01);
#endif    

    using CT = typename bg::coordinate_type<Point>::type;
    geo_buffer_accurate_area<Formula, Spheroid, CT> strategy(spheroid);

    bg::strategy::buffer::geographic_point_circle<Formula, Spheroid> circle(spheroid, points_per_circle);

    // All are ignored for points
    bg::strategy::buffer::geographic_side_straight<Formula, Spheroid> side(spheroid);
    bg::strategy::buffer::geographic_join_miter<Formula, Spheroid> join_miter(spheroid);
    bg::strategy::buffer::geographic_join_round<Formula, Spheroid> join_round(spheroid);
    bg::strategy::buffer::geographic_end_round<Formula, Spheroid> end_round(spheroid);

    test_one_geo<Point, polygon>(label, torg, strategy, side, circle, join_round, end_round, expected_area, 100.0, settings);
}

template <typename test_type>
void test_all()
{
    using point_t = bg::model::point<test_type, 2, bg::cs::geographic<bg::degree>>;
    using strategy = bg::strategy::andoyer;

    // Use the default spheroid
    bg::srs::spheroid<test_type> def_spheroid;
    test_linestring<strategy, true, point_t>("line_def", def_spheroid, 8046.26, 8143.37);
    test_point<strategy, true, point_t>("point_def", def_spheroid, 31414.36);

    // Call it with a quite different spheroid (a near sphere), this changes internal geographic calculations
    // and should result in different areas. Using CSV creation, it's visible in QGis.
    bg::srs::spheroid<test_type> alt_spheroid(6378000.0, 6375000.0);
    test_linestring<strategy, true, point_t>("line_alt", alt_spheroid, 8030.53, 8127.61);
    test_point<strategy, true, point_t>("point_alt", alt_spheroid, 31414.33);
}

int test_main(int, char* [])
{
    BoostGeometryWriteTestConfiguration();

    // There are several issues with float type such as invalid geometries as output and assertion fails
    //test_all<float>();
    test_all<default_test_type>();

    return 0;
}
