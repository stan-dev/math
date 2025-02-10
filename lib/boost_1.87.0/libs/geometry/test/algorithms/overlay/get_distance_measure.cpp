// Boost.Geometry
// Unit Test

// Copyright (c) 2023 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <geometry_test_common.hpp>

#include <boost/geometry/algorithms/detail/overlay/get_distance_measure.hpp>
#include <boost/geometry/strategy/cartesian/side_robust.hpp>
#include <boost/geometry/strategies/strategies.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/io/wkt/wkt.hpp>

// #define BOOST_GEOMETRY_TEST_WITH_COUT
// #define BOOST_GEOMETRY_TEST_FAILURES

template <typename Point>
void do_test(std::string const& case_id,
             Point const& s1,
             Point const& s2,
             Point const& p,
             int expected_side,
             bool ignore_failure = false)
{
    using coor_t = typename bg::coordinate_type<Point>::type;
    typename bg::strategies::relate::services::default_strategy
        <
            Point, Point
        >::type strategy;

    auto const dm = bg::detail::get_distance_measure(s1, s2, p, strategy);
    auto const dm_side = dm.measure < 0.0 ? -1 : dm.measure > 0.0 ? 1 : 0;
    auto const tr_side = strategy.side().apply(s1, s2, p);
    auto const rob_side = bg::strategy::side::side_robust<coor_t>::apply(s1, s2, p);

    #if defined(BOOST_GEOMETRY_TEST_WITH_COUT)
    std::cout << " " << case_id << " " << string_from_type<coor_t>::name()
        << std::setprecision(20)
        << " " << dm.measure << " " << dm_side << " " << tr_side << " " << rob_side
        << (dm_side != rob_side ? " [*** DM WRONG]" : "")
        << (tr_side != rob_side ? " [*** TR WRONG]" : "")
        << std::endl;
    #endif

    BOOST_CHECK_MESSAGE(expected_side == -9 || expected_side == dm_side,
                        "Case: " << case_id
                        << " ctype: " << string_from_type<coor_t>::name()
                        << " expected: " << expected_side
                        << " detected: " << dm_side);

    // This is often wrong for float, and sometimes for double
    BOOST_CHECK_MESSAGE(ignore_failure || tr_side == dm_side,
                        "Case: " << case_id
                        << " ctype: " << string_from_type<coor_t>::name()
                        << " tr_side: " << tr_side
                        << " dm_side: " << dm_side);

    // This is always guaranteed for the tested values
    BOOST_CHECK_MESSAGE(tr_side == rob_side,
                        "Case: " << case_id
                        << " ctype: " << string_from_type<coor_t>::name()
                        << " tr_side: " << tr_side
                        << " rob_side: " << rob_side);
}

template <typename Point>
void test_get_distance_measure()
{
    using coor_t = typename bg::coordinate_type<Point>::type;

    do_test<Point>("simplex_coll", {1.0, 0.0}, {1.0, 1.0}, {1.0, 0.5}, 0);
    do_test<Point>("simplex_left", {1.0, 0.0}, {1.0, 1.0}, {0.9, 0.5}, 1);
    do_test<Point>("simplex_right", {1.0, 0.0}, {1.0, 1.0}, {1.1, 0.5}, -1);

    bool const is_float = std::is_floating_point<coor_t>::value && sizeof(coor_t) == 4;
    bool const is_double = std::is_floating_point<coor_t>::value && sizeof(coor_t) == 8;

    // The issue 1183 where get_distance_measure failed for these coordinates.
    std::string const case_id = "issue_1183_";
    Point const p1{38902.349206128216, 6721371.1493254723};
    Point const p2{38937.993505971914, 6721407.9151819283};
    Point const q1 = p1;
    Point const q2{38960.647313876834, 6721431.2817974398};
    {
        do_test<Point>(case_id + "p", p1, p2, q2, 1, is_float);
        do_test<Point>(case_id + "q", q1, q2, p1, 0);
    }

    double const test_epsilon = 1.0e-9;

    // Walk along x axis to get the switch from left to right (between 2 and 3)
    for (int i = -5; i <= 15; i++)
    {
#if defined(BOOST_GEOMETRY_TEST_FAILURES)
        bool const ignore_failure = false;
#else
        bool const ignore_failure = is_float || (is_double && i >= 3 && i <= 12);
#endif
        double const v = i / 10.0;
        Point q2a = q2;
        bg::set<0>(q2a, bg::get<0>(q2) + v * test_epsilon);
        do_test<Point>(case_id + std::to_string(i), p1, p2, q2a, -9, ignore_failure);
    }
    // Focus on the switch from left to right (between 0.251 and 0.252)
    for (int i = 250; i <= 260; i++)
    {
        double const v = i / 1000.0;
        Point q2a = q2;
        bg::set<0>(q2a, bg::get<0>(q2) + v * test_epsilon);
        do_test<Point>(case_id + std::to_string(i), p1, p2, q2, -9, is_float || is_double);
    }
}

int test_main(int, char* [])
{
    using fp = bg::model::point<float, 2, bg::cs::cartesian>;
    using dp = bg::model::point<double, 2, bg::cs::cartesian>;
    using ep = bg::model::point<long double, 2, bg::cs::cartesian>;

    test_get_distance_measure<fp>();
    test_get_distance_measure<dp>();
    test_get_distance_measure<ep>();

    return 0;
}
