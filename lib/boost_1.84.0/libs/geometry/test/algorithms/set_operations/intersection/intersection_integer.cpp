// Boost.Geometry
// Unit Test

// Copyright (c) 2007-2023 Barend Gehrels, Amsterdam, the Netherlands.

// This file was modified by Oracle on 2015-2022.
// Modifications copyright (c) 2015-2022, Oracle and/or its affiliates.
// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "test_intersection.hpp"
#include <algorithms/test_overlay.hpp>

#include <algorithms/overlay/overlay_cases.hpp>

#define TEST_INTERSECTION(caseid, clips, points, area) \
    (test_one<Polygon, Polygon, Polygon>) \
    ( #caseid, caseid[0], caseid[1], clips, points, area, settings)

namespace
{
    std::string rectangles[2] =
    {
        "POLYGON((1000 1000,2000 1000,2000 2000,1000 2000))",
        "POLYGON((500 500,1500 500,1500 1500,500 1500))"
    };

    std::string start_within[2] =
    {
        "POLYGON((1000 1000,2000 1000,2000 2000,1000 2000))",
        "POLYGON((1250 1250,1500 750,1750 1250,1750 1750))"
    };
}


template <typename Polygon>
void test_areal()
{
    static bool is_open = bg::closure<Polygon>::value == bg::open;

    ut_settings settings;
    settings.test_point_count = true;

    TEST_INTERSECTION(rectangles, 1, is_open ? 4 : 5, 250000.0);
    TEST_INTERSECTION(start_within, 1, is_open ? 5 : 6, 218750.0);
    TEST_INTERSECTION(issue_1184, 1, is_open ? 3 : 4, 156.0);
}

template <typename P>
void test_all()
{
    using polygon = bg::model::polygon<P>;
    using polygon_ccw = bg::model::polygon<P, false>;
    using polygon_open = bg::model::polygon<P, true, false>;
    using polygon_ccw_open = bg::model::polygon<P, false, false>;
    test_areal<polygon>();
    test_areal<polygon_ccw>();
    test_areal<polygon_open>();
    test_areal<polygon_ccw_open>();
}

template <typename CoordinateType>
void test_ticket_10868(/*std::string const& wkt_out*/)
{
    using point_type = bg::model::point<CoordinateType, 2, bg::cs::cartesian>;
    using polygon_type = bg::model::polygon
        <
            point_type, /*ClockWise*/false, /*Closed*/false
        >;
    using multipolygon_type = bg::model::multi_polygon<polygon_type>;

    polygon_type polygon1;
    polygon_type polygon2;
    bg::read_wkt(ticket_10868[0], polygon1);
    bg::read_wkt(ticket_10868[1], polygon2);

    multipolygon_type multipolygon_out;
    bg::intersection(polygon1, polygon2, multipolygon_out);
    std::stringstream stream;
    stream << bg::wkt(multipolygon_out);

    // BOOST_CHECK_EQUAL(stream.str(), wkt_out);

    test_one<polygon_type, polygon_type, polygon_type>("ticket_10868",
        ticket_10868[0], ticket_10868[1],
        1, 7, 20266195244586.0);
}



int test_main(int, char* [])
{
    BoostGeometryWriteTestConfiguration();

    test_all<bg::model::d2::point_xy<std::int32_t> >();


#if defined(BOOST_GEOMETRY_TEST_FAILURES)
    // ticket #10868 still fails for 32-bit integers
    test_ticket_10868<std::int32_t>("MULTIPOLYGON(((33520458 6878575,33480192 14931538,31446819 18947953,30772384 19615678,30101303 19612322,30114725 16928001,33520458 6878575)))");

    test_ticket_10868<std::int64_t>("MULTIPOLYGON(((33520458 6878575,33480192 14931538,31446819 18947953,30772384 19615678,30101303 19612322,30114725 16928001,33520458 6878575)))");
#endif

 #if defined(BOOST_GEOMETRY_TEST_FAILURES)
    // Ticket #10868 was normally not run for other types
    // It results in a different area (might be fine)
    // and it reports self intersections.
    if (BOOST_GEOMETRY_CONDITION(sizeof(long) * CHAR_BIT >= 64))
    {
        test_ticket_10868<long>("MULTIPOLYGON(((33520458 6878575,33480192 14931538,31446819 18947953,30772384 19615678,30101303 19612322,30114725 16928001,33520458 6878575)))");
    }

    test_ticket_10868<long long>("MULTIPOLYGON(((33520458 6878575,33480192 14931538,31446819 18947953,30772384 19615678,30101303 19612322,30114725 16928001,33520458 6878575)))");
 #endif

    return 0;
}
