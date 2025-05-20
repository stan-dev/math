// Boost.Geometry (aka GGL, Generic Geometry Library)
// Unit Test

// Copyright (c) 2007-2024 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#include <geometry_test_common.hpp>

#include <boost/geometry/algorithms/area.hpp>
#include <boost/geometry/algorithms/assign.hpp>
#include <boost/geometry/algorithms/make.hpp>
#include <boost/geometry/algorithms/expand.hpp>
#include <boost/geometry/algorithms/within.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/io/wkt/wkt.hpp>
#include <boost/geometry/util/numeric_cast.hpp>
#include <boost/geometry/util/rational.hpp>

void test_coordinate_cast(std::string const& s, int expected_nom, int expected_denom)
{
    boost::rational<int> const a = bg::detail::coordinate_cast<boost::rational<int> >::apply(s);
    BOOST_CHECK_EQUAL(a.numerator(), expected_nom);
    BOOST_CHECK_EQUAL(a.denominator(), expected_denom);
}

void test_numeric_cast()
{
    boost::rational<int> const r1(3, 4);
    BOOST_CHECK_CLOSE(bg::util::numeric_cast<double>(r1), 0.75, 0.00001);

    boost::rational<int> const r2(10, 4);
    BOOST_CHECK_CLOSE(bg::util::numeric_cast<double>(r2), 2.5, 0.00001);
    BOOST_CHECK_EQUAL(bg::util::numeric_cast<int>(r2), 2);
}

template <typename T>
void test_bounds()
{
    using coordinate_t = boost::rational<T>;
    using point_t = bg::model::point<coordinate_t, 2, bg::cs::cartesian>;

    auto const lowest = bg::util::bounds<coordinate_t>::lowest();
    auto const highest = bg::util::bounds<coordinate_t>::highest();

    BOOST_CHECK_MESSAGE(lowest < highest,
        "Lowest should be smaller than highest, lowest: " << lowest << " highest: " << highest);
}

// Tests box-related functionality, which depends on geometry::util::bounds
// specialization for Boost.Rational
template <typename T>
void test_box()
{
    using coordinate_t = boost::rational<T>;
    using point_t = bg::model::point<coordinate_t, 2, bg::cs::cartesian>;
    using box_t = bg::model::box<point_t>;

    box_t box;
    bg::assign_inverse(box);

    point_t south_west, north_east;
    bg::detail::assign_point_from_index<0>(box, south_west);
    bg::detail::assign_point_from_index<1>(box, north_east);

    BOOST_CHECK_MESSAGE(bg::get<0>(south_west) > bg::get<0>(north_east),
        "Bounding box should be inversed. Now x-min: " << bg::get<0>(south_west)
            << " x-max: " << bg::get<0>(north_east)
            << " " << bg::wkt(box));

    BOOST_CHECK_MESSAGE(bg::get<1>(south_west) > bg::get<1>(north_east),
        "Bounding box should be inversed. Now y-min: " << bg::get<1>(south_west)
            << " y-max: " << bg::get<1>(north_east)
            << " " << bg::wkt(box));

    // Test specifically for points larger than 0, because without specialization Boost.Rational
    // will return (0,1) (== 0) by default and code will compile but give wrong results.
    bg::expand(box, bg::make<point_t>(4, 4));
    bg::expand(box, bg::make<point_t>(8, 8));

    // Test within (without specialization, both points are within the box)
    auto const point1 = bg::make<point_t>(6, 6);
    auto const point2 = bg::make<point_t>(2, 2);
    BOOST_CHECK_MESSAGE(bg::within(point1, box),
        "Point " << bg::wkt(point1) << " is not within the box " << bg::wkt(box));
    BOOST_CHECK_MESSAGE(! bg::within(point2, box),
        "Point " << bg::wkt(point2) << " is within the box " << bg::wkt(box));

    // Test area (without specialization, it will be 64)
    auto const area = bg::util::numeric_cast<T>(bg::area(box));
    T const expected_area = 16;
    BOOST_CHECK_EQUAL(expected_area, area);
}

void test_select_most_precise()
{
    using rational1_t = boost::rational<std::int32_t>;
    using rational2_t = boost::rational<std::int64_t>;

    using t1 = bg::select_most_precise<double, rational1_t>::type;
    using t2 = bg::select_most_precise<double, rational2_t>::type;
    using t12 = bg::select_most_precise<rational1_t, rational2_t>::type;

    BOOST_CHECK((std::is_same<t1, rational1_t>::value));
    BOOST_CHECK((std::is_same<t2, rational2_t>::value));
    BOOST_CHECK((std::is_same<t12, rational2_t>::value));
}

void test_wkt(std::string const& wkt, std::string const expected_wkt)
{
    bg::model::point<boost::rational<int>, 2, bg::cs::cartesian> p;
    bg::read_wkt(wkt, p);
    std::ostringstream out;
    out << bg::wkt(p);

    BOOST_CHECK_EQUAL(out.str(), expected_wkt);
}

int test_main(int, char* [])
{
    test_coordinate_cast("0", 0, 1);
    test_coordinate_cast("1", 1, 1);
    test_coordinate_cast("-1", -1, 1);
    test_coordinate_cast("-0.5", -1, 2);
    test_coordinate_cast("-1.5", -3, 2);
    test_coordinate_cast("0.5", 1, 2);
    test_coordinate_cast("1.5", 3, 2);
    test_coordinate_cast("2.12345", 42469, 20000);
    test_coordinate_cast("1.", 1, 1);

    test_coordinate_cast("3/2", 3, 2);
    test_coordinate_cast("-3/2", -3, 2);

    test_numeric_cast();

    test_bounds<std::int16_t>();
    test_bounds<std::int32_t>();
    test_bounds<std::int64_t>();

    test_box<std::int64_t>();

    test_select_most_precise();

    test_wkt("POINT(1.5 2.75)", "POINT(3/2 11/4)");
    test_wkt("POINT(3/2 11/4)", "POINT(3/2 11/4)");
    test_wkt("POINT(-1.5 2.75)", "POINT(-3/2 11/4)");
    test_wkt("POINT(-3/2 11/4)", "POINT(-3/2 11/4)");

    return 0;
}
