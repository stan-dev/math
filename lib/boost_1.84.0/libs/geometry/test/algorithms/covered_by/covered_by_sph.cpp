// Boost.Geometry

// Copyright (c) 2016 Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "test_covered_by.hpp"

#include <algorithms/overlay/overlay_cases.hpp>
#include <algorithms/overlay/multi_overlay_cases.hpp>

#include <boost/geometry/geometries/geometries.hpp>

template <typename P>
void test_polygon_polygon()
{
    test_geometry<ring<P>, ring<P>>(case_1[0], case_1[1], false);
    test_geometry<ring<P>, poly<P>>(case_1[0], case_1[1], false);

    test_geometry<poly<P>, poly<P>>(case_1[0], case_1[1], false);
    test_geometry<poly<P>, poly<P>>(case_2[0], case_2[1], false);
    test_geometry<poly<P>, poly<P>>(case_3_sph[0], case_3_sph[1], true);
    test_geometry<poly<P>, poly<P>>(case_3_2_sph[0], case_3_2_sph[1], true);
    test_geometry<poly<P>, poly<P>>(case_4[0], case_4[1], false);
    test_geometry<poly<P>, poly<P>>(case_5[0], case_5[1], false);
    test_geometry<poly<P>, poly<P>>(case_6_sph[0], case_6_sph[1], false);
    test_geometry<poly<P>, poly<P>>(case_6_sph[1], case_6_sph[0], true);

    test_geometry<poly<P>, poly<P>>(case_7[0], case_7[1], false);
    test_geometry<poly<P>, poly<P>>(case_8_sph[0], case_8_sph[1], false);
    test_geometry<poly<P>, poly<P>>(case_9_sph[0], case_9_sph[1], false);
    test_geometry<poly<P>, poly<P>>(case_10_sph[0], case_10_sph[1], false);
    test_geometry<poly<P>, poly<P>>(case_11_sph[0], case_11_sph[1], false);
    test_geometry<poly<P>, poly<P>>(case_11_sph[1], case_11_sph[0], true);
    test_geometry<poly<P>, poly<P>>(case_12[0], case_12[1], false);

    test_geometry<poly<P>, poly<P>>(case_13_sph[0], case_13_sph[1], false);
    test_geometry<poly<P>, poly<P>>(case_14_sph[0], case_14_sph[1], false);
    test_geometry<poly<P>, poly<P>>(case_15_sph[0], case_15_sph[1], false);
    test_geometry<poly<P>, poly<P>>(case_16_sph[0], case_16_sph[1], false);
    test_geometry<poly<P>, poly<P>>(case_17_sph[0], case_17_sph[1], false);
    test_geometry<poly<P>, poly<P>>(case_17_sph[1], case_17_sph[0], true);
    test_geometry<poly<P>, poly<P>>(case_18_sph[0], case_18_sph[1], false);
    test_geometry<poly<P>, poly<P>>(case_18_sph[1], case_18_sph[0], true);
}

template <typename P>
void test_polygon_multi_polygon()
{
    test_geometry<ring<P>, mpoly<P>>(case_1[0], case_multi_2[0], false);
    test_geometry<poly<P>, mpoly<P>>(case_2[0], case_multi_2[0], false);
}

template <typename P>
void test_multi_polygon_multi_polygon()
{
    test_geometry<mpoly<P>, mpoly<P>>(case_multi_2[0], case_multi_2[1], false);
}

template <typename P>
void test_linestring_polygon()
{
    test_geometry<ls<P>, poly<P>>("LINESTRING(11 0,11 10)", "POLYGON((0 0,0 10,10 10,10 0,0 0))", false);
    test_geometry<ls<P>, ring<P>>("LINESTRING(11 0,11 10)", "POLYGON((0 0,0 10,10 10,10 0,0 0))", false);
    test_geometry<ls<P>, poly<P>>("LINESTRING(0 0,10 10)", "POLYGON((0 0,0 10,10 10,10 0,0 0))", true);
    test_geometry<ls<P>, poly<P>>("LINESTRING(5 0,5 5,10 5)", "POLYGON((0 0,0 10,10 10,10 0,0 0))", true);
    test_geometry<ls<P>, poly<P>>("LINESTRING(5 1,5 5,9 5)", "POLYGON((0 0,0 10,10 10,10 0,0 0))", true);
    test_geometry<ls<P>, poly<P>>("LINESTRING(11 1,11 5)", "POLYGON((0 0,0 10,10 10,10 0,0 0))", false);

    test_geometry<ls<P>, poly<P>>("LINESTRING(9 1,10 5,9 9)",
                                  "POLYGON((0 0,0 10,10 10,10 0,0 0),(10 5,2 8,2 2,10 5))",
                                  true);

    test_geometry<ls<P>, poly<P>>("LINESTRING(9 1,10 5,9 9,1 9,1 1,9 1)",
                                  "POLYGON((0 0,0 10,10 10,10 0,0 0),(10 5,2 8,2 2,10 5))",
                                  true);

    test_geometry<ls<P>, poly<P>>("LINESTRING(0 0,10 0,10 10,0 10,0 0)",
                                  "POLYGON((0 0,0 10,10 10,10 0,0 0))",
                                  true);
}

template <typename P>
void test_linestring_multi_polygon()
{
    test_geometry<ls<P>, mpoly<P>>("LINESTRING(10 1,10 5,10 9)",
                                   "MULTIPOLYGON(((0 20,0 30,10 30,10 20,0 20)),((0 0,0 10,10 10,10 0,0 0),(10 5,2 8,2 2,10 5)))",
                                   true);
}

template <typename P>
void test_multi_linestring_polygon()
{
    test_geometry<mls<P>, poly<P>>("MULTILINESTRING((11 11, 20 20),(5 7, 4 1))",
                                   "POLYGON((0 0,0 10,10 10,10 0,0 0),(2 2,4 2,4 4,2 4,2 2))",
                                   false);

    test_geometry<mls<P>, ring<P>>("MULTILINESTRING((6 6,15 15),(0 0, 7 7))",
                                   "POLYGON((5 5,5 15,15 15,15 5,5 5))",
                                   false);

    test_geometry<mls<P>, poly<P>>("MULTILINESTRING((3 10.031432746397092, 1 5, 1 10.013467818052765, 3 4, 7 8, 6 10.035925377760330, 10 2))",
                                   "POLYGON((0 0,0 10,10 10,10 0,0 0))",
                                   true);
}

template <typename P>
void test_multi_linestring_multi_polygon()
{
    test_geometry<mls<P>, mpoly<P>>("MULTILINESTRING((0 0,10 0,10 10,0 10,0 0),(2 2,5 5,2 8,2 2))",
                                    "MULTIPOLYGON(((0 0,0 10,10 10,10 0,0 0),(2 2,5 5,2 8,2 2)))",
                                    true);

    test_geometry<mls<P>, mpoly<P>>("MULTILINESTRING((0 0,10 0,10 10),(10 10,0 10,0 0),(20 20,50 50,20 80,20 20))",
                                    "MULTIPOLYGON(((0 0,0 10,10 10,10 0,0 0)))",
                                    false);

    test_geometry<mls<P>, mpoly<P>>("MULTILINESTRING((5 -2,4 -2,5 0),(5 -2,6 -2,5 0))",
                                    "MULTIPOLYGON(((5 0,0 5,10 5,5 0)),((5 0,10 -5,0 -5,5 0)))",
                                    true);
}

template <typename P>
void test_linestring_linestring()
{
    test_geometry<ls<P>, ls<P>>("LINESTRING(0 0, 2 2, 3 2)", "LINESTRING(0 0, 2 2, 3 2)", true);

    test_geometry<ls<P>, ls<P>>("LINESTRING(1 0,2 2,2 3)", "LINESTRING(0 0, 2 2, 3 2)", false);
}

template <typename P>
void test_linestring_multi_linestring()
{
    test_geometry<ls<P>, mls<P>>("LINESTRING(0 0,10 0)",
                                "MULTILINESTRING((1 0,2 0),(1 1,2 1))",
                                false);

    test_geometry<ls<P>, mls<P>>("LINESTRING(0 0,5 0,5 5,0 5,0 0)",
                                 "MULTILINESTRING((5 5,0 5,0 0),(0 0,5 0,5 5))",
                                 true);
}

template <typename P>
void test_multi_linestring_multi_linestring()
{
    test_geometry<mls<P>, mls<P>>("MULTILINESTRING((0 0,0 0,18 0,18 0,19 0,19 0,19 0,30 0,30 0))",
                                  "MULTILINESTRING((0 10,5 0,20 0,20 0,30 0))",
                                  false);
}

template <typename P>
void test_point_polygon()
{
    // https://svn.boost.org/trac/boost/ticket/9162
    test_geometry<P, poly<P>>("POINT(0 90)",
                              "POLYGON((0 80,-90 80, -180 80, 90 80, 0 80))",
                              true);
    test_geometry<P, poly<P>>("POINT(-120 21)",
                              "POLYGON((30 0,30 30,90 30, 90 0, 30 0))",
                              false);
    // extended
    test_geometry<P, poly<P>>("POINT(0 90)",
                              "POLYGON((0 80, 0 81, -90 80, -180 80, 90 80, 0 80))",
                              true);
    test_geometry<P, poly<P>>("POINT(0 90)",
                              "POLYGON((0 80, -90 80, -90 81, -180 80, 90 80, 0 80))",
                              true);
    test_geometry<P, poly<P>>("POINT(0 89)",
                              "POLYGON((0 80,-90 80, -180 80, 90 80, 0 80))",
                              true);
    test_geometry<P, poly<P>>("POINT(-180 89)",
                              "POLYGON((0 80,-90 80, -180 80, 90 80, 0 80))",
                              true);
    test_geometry<P, poly<P>>("POINT(0 -90)",
                              "POLYGON((0 -80,90 -80, -180 -80, -90 -80, 0 -80))",
                              true);
    test_geometry<P, poly<P>>("POINT(0 -89)",
                              "POLYGON((0 -80,90 -80, -180 -80, -90 -80, 0 -80))",
                              true);
    test_geometry<P, poly<P>>("POINT(1 -90)",
                              "POLYGON((0 -80,90 -80, -180 -80, -90 -80, 0 -80))",
                              true);
    test_geometry<P, poly<P>>("POINT(1 -89)",
                              "POLYGON((0 -80,90 -80, -180 -80, -90 -80, 0 -80))",
                              true);
    test_geometry<P, poly<P>>("POINT(1 90)",
                              "POLYGON((0 80,-90 80, -180 80, 90 80, 0 80))",
                              true);
    test_geometry<P, poly<P>>("POINT(1 90)",
                              "POLYGON((0 80,-90 80, -180 80, 90 80, 0 80))",
                              true);



    // MySQL report 08.2017
    test_geometry<P, poly<P>>("POINT(-179 0)",
                              "POLYGON((0 0, 0 2, 2 0, 0 -2, 0 0))",
                              false);
    // extended
    test_geometry<P, poly<P>>("POINT(179 0)",
                              "POLYGON((0 0, 0 2, 2 0, 0 -2, 0 0))",
                              false);
    test_geometry<P, poly<P>>("POINT(180 0)",
                              "POLYGON((0 0, 0 2, 2 0, 0 -2, 0 0))",
                              false);
    test_geometry<P, poly<P>>("POINT(-179 0)",
                              "POLYGON((-10 -10, -10 10, 10 10, 10 -10, -10 10))",
                              false);
    test_geometry<P, poly<P>>("POINT(179 0)",
                              "POLYGON((-10 -10, -10 10, 10 10, 10 -10, -10 10))",
                              false);
    test_geometry<P, poly<P>>("POINT(-179 0)",
                              "POLYGON((0 0, 0 1, 1 0, 0 -1, 0 0))",
                              false);
    test_geometry<P, poly<P>>("POINT(179 0)",
                              "POLYGON((0 0, 0 1, 1 0, 0 -1, 0 0))",
                              false);
}


template <typename P>
void test_all()
{
    test_polygon_polygon<P>();
    test_polygon_multi_polygon<P>();
    test_multi_polygon_multi_polygon<P>();

    test_linestring_polygon<P>();
    test_linestring_multi_polygon<P>();
    test_multi_linestring_polygon<P>();
    test_multi_linestring_multi_polygon<P>();

    test_linestring_linestring<P>();
    test_linestring_multi_linestring<P>();
    test_multi_linestring_multi_linestring<P>();

    test_point_polygon<P>();
}


int test_main( int , char* [] )
{
    test_all<bg::model::point<double, 2, bg::cs::spherical_equatorial<bg::degree> > >();

    return 0;
}
