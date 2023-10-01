// Boost.Geometry
// Unit Test

// Copyright (c) 2022 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2023 Adam Wulkiewicz, Lodz, Poland.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <geometry_test_common.hpp>

#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/strategies/geographic/buffer_join_round.hpp>
#include <boost/geometry/io/wkt/wkt.hpp>
#include <boost/geometry/algorithms/area.hpp>
#include <boost/geometry/algorithms/distance.hpp>
#include <boost/geometry/algorithms/equals.hpp>


template<typename Point, typename T>
inline Point generate_point(Point const& point, T const& distance, T const& azimuth)
{
    constexpr bool enable_coordinates = true;

    using calc_t = typename bg::coordinate_type<Point>::type;
    using direct_t = bg::strategy::andoyer::direct
        <
            calc_t, enable_coordinates, false, false, false
        >;

    bg::srs::spheroid<calc_t> spheroid;
    auto const angle_rad = bg::math::wrap_azimuth_in_radian(azimuth * bg::math::d2r<calc_t>());
    auto const d = direct_t::apply(bg::get_as_radian<0>(point), bg::get_as_radian<1>(point),
                                    distance, angle_rad, spheroid);
    Point p;
    bg::set_from_radian<0>(p, d.lon2);
    bg::set_from_radian<1>(p, d.lat2);
    return p;
}

template <typename JoinStrategy, typename P, typename T>
void test_join(std::string const& case_id, JoinStrategy const& join, P const& vertex,
               T const& buffer_distance, T const& angle1, T const& angle2,
               std::size_t expected_size)
{
    boost::ignore_unused(case_id);

    // Use a deque to be able to use push_front
    bg::model::ring<P, true, true, std::deque> output_ring;

    // IP is not used for geographic join
    P ip{0.0, 0.0};

    auto const perp1 = generate_point(vertex, buffer_distance, angle1);
    auto const perp2 = generate_point(vertex, buffer_distance, angle2);

    join.apply(ip, vertex, perp1, perp2, buffer_distance, output_ring);

    BOOST_CHECK_EQUAL(expected_size , output_ring.size());

    // All the generated points should be located
    // at or close to the specified buffer distance from the vertex
    for (const auto& p : output_ring)
    {
        auto const d = bg::distance(vertex, p);
        auto const fraction = d / buffer_distance;
        BOOST_CHECK_MESSAGE(fraction > 0.99 && fraction < 1.01,
                "Unexpected distance = " << d << " fraction=" << fraction);
    }

#ifdef TEST_WITH_CSV
    // Optionally create output for QGis to inspect the result visually
    std::string path = "/tmp/csv/";
    std::ofstream out(path + case_id + ".csv");
    if (out.good())
    {
        out << std::setprecision(16) << "id;role;wkt" << std::endl;
        int id = 1;
        out << id++ << ";1;" << bg::wkt(vertex) << std::endl;
        out << id++ << ";2;" << bg::wkt(perp1) << std::endl;
        out << id++ << ";2;" << bg::wkt(perp2) << std::endl;
        for (const auto& p : output_ring)
        {
            out << id++ << ";3;" << bg::wkt(p) << std::endl;
        }
    }
#endif
}

template <typename P>
void test_all()
{
    // Test point on the Dam of Amsterdam.
    P const dam{4.8924323, 52.3731145};

    bg::strategy::buffer::geographic_join_round<> join(72);
    test_join("part1", join, dam, 50.0, 45.0, 135.0, 20);
    test_join("part2", join, dam, 25.0, 135.0, 310.0, 36);
    test_join("part3", join, dam, 40.0, 310.0, 45.0, 20);

    bg::strategy::buffer::geographic_join_round<> join360(360);
    test_join("large_part1", join360, dam, 10000.0, 0.0, 90.0, 92);
    test_join("large_part2", join360, dam, 10100.0, 90.0, 269.0, 180);
    test_join("large_part3", join360, dam, 10200.0, 269.0, 10.0, 102);
}

int test_main(int, char *[])
{
    test_all<bg::model::point<double, 2, bg::cs::geographic<bg::degree>>>();
    return 0;
}
