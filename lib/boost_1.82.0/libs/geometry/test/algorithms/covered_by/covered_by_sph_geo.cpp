// Boost.Geometry

// Copyright (c) 2016-2022 Oracle and/or its affiliates.

// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "test_covered_by.hpp"


#include <boost/geometry/geometries/geometries.hpp>


template <typename P>
void test_point_box()
{
    test_geometry<P, box<P>>("POINT(0 0)",    "BOX(0 0, 1 1)", true);
    test_geometry<P, box<P>>("POINT(1 1)",    "BOX(0 0, 2 2)", true);

    test_geometry<P, box<P>>("POINT(180 1)",  "BOX(170 0, 190 2)", true);
    test_geometry<P, box<P>>("POINT(-180 1)", "BOX(170 0, 190 2)", true);
    test_geometry<P, box<P>>("POINT(180 1)",  "BOX(170 0, 180 2)", true);
    test_geometry<P, box<P>>("POINT(-180 1)", "BOX(170 0, 180 2)", true);
    test_geometry<P, box<P>>("POINT(179 1)",  "BOX(170 0, 190 2)", true);
    test_geometry<P, box<P>>("POINT(-179 1)", "BOX(170 0, 190 2)", true);
    test_geometry<P, box<P>>("POINT(179 1)",  "BOX(170 0, 180 2)", true);
    test_geometry<P, box<P>>("POINT(-179 1)", "BOX(170 0, 180 2)", false);
    test_geometry<P, box<P>>("POINT(169 1)", "BOX(170 0, 180 2)", false);

    // https://svn.boost.org/trac/boost/ticket/12412
    test_geometry<P, box<P>>("POINT(-0.127592 51.7)", "BOX(-2.08882 51.5034, -0.127592 51.9074)", true);
    // and related
    test_geometry<P, box<P>>("POINT(-2.08882 51.7)", "BOX(-2.08882 51.5034, -0.127592 51.9074)", true);
    test_geometry<P, box<P>>("POINT(0.127592 51.7)", "BOX(0.127592 51.5034, 2.08882 51.9074)", true);
    test_geometry<P, box<P>>("POINT(2.08882 51.7)", "BOX(0.127592 51.5034, 2.08882 51.9074)", true);

    test_geometry<P, box<P>>("POINT(179.08882 1)", "BOX(179.08882 0, 538.127592 2)", true);
    test_geometry<P, box<P>>("POINT(178.127592 1)", "BOX(179.08882 0, 538.127592 2)", true);
    test_geometry<P, box<P>>("POINT(179.08882 1)", "BOX(179.08882 0, 182.127592 2)", true);
    test_geometry<P, box<P>>("POINT(-177.872408 1)", "BOX(179.08882 0, 182.127592 2)", true);
}

template <typename P>
void test_multi_point_box()
{
    test_geometry<mpt<P>, box<P>>("MULTIPOINT(0 0,1 1)",    "BOX(0 0, 1 1)", true);
    test_geometry<mpt<P>, box<P>>("MULTIPOINT(1 1,3 3)",    "BOX(0 0, 2 2)", false);
}

template <typename P>
void test_areal_box()
{
    test_geometry<ring<P>, box<P>>("POLYGON((0 0,0 3,3 3,3 0,0 0))", "BOX(0 0,4 4)", true);
    test_geometry<ring<P>, box<P>>("POLYGON((0 0,0 3,3 3,5 0,0 0))", "BOX(0 0,4 4)", false);
    test_geometry<poly<P>, box<P>>("POLYGON((0 0,0 3,3 3,3 0,0 0))", "BOX(0 0,4 4)", true);
    test_geometry<poly<P>, box<P>>("POLYGON((0 0,0 3,3 3,5 0,0 0))", "BOX(0 0,4 4)", false);
    // the following is true for cartesian but not for spherical or geographic
    // because in non-cartesian boxes the horizontal edges are not geodesics
    test_geometry<mpoly<P>, box<P>>("MULTIPOLYGON(((0 0,0 3,3 3,3 0,0 0)),((4 4,4 7,7 7,4 7,4 4)))",
                                  "BOX(0 0,7 7)", false);
    test_geometry<mpoly<P>, box<P>>("MULTIPOLYGON(((0 0,0 3,3 3,3 0,0 0)),((4 4,4 6.5,6.5 6.5,4 6.5,4 4)))",
                                  "BOX(0 0,7 7)", true);
    test_geometry<mpoly<P>, box<P>>("MULTIPOLYGON(((0 0,0 3,3 3,5 0,0 0)),((4 4,4 7,7 7,4 7,4 4)))",
                                  "BOX(0 0,4 4)", false);

}

template <typename P>
void test_box_box()
{
    test_geometry<box<P>, box<P>>("BOX(0 0, 1 1)", "BOX(0 0, 1 1)", true);

    test_geometry<box<P>, box<P>>("BOX(-170 0,-160 1)", "BOX(-180 0, 180 1)", true);
    test_geometry<box<P>, box<P>>("BOX(-170 0,-160 1)", "BOX(170 0, 200 1)",  true);
    test_geometry<box<P>, box<P>>("BOX(-170 0,-150 1)", "BOX(170 0, 200 1)",  false);
    test_geometry<box<P>, box<P>>("BOX(0 0,1 1)",       "BOX(170 0, 370 1)",  true);
    test_geometry<box<P>, box<P>>("BOX(0 0,10 1)",      "BOX(170 0, 370 1)",  true);
    test_geometry<box<P>, box<P>>("BOX(-180 0,10 1)",   "BOX(170 0, 370 1)",  true);
    test_geometry<box<P>, box<P>>("BOX(-180 0,20 1)",   "BOX(170 0, 370 1)",  false);
    test_geometry<box<P>, box<P>>("BOX(10 0,20 1)",     "BOX(170 0, 370 1)",  false);
    test_geometry<box<P>, box<P>>("BOX(160 0,180 1)",   "BOX(170 0, 370 1)",  false);

    test_geometry<box<P>, box<P>>("BOX(-180 0,-170 1)", "BOX(180 0, 190 1)",  true); // invalid?
    test_geometry<box<P>, box<P>>("BOX(-180 0,-170 1)", "BOX(180 0, 191 1)",  true); // invalid?
    test_geometry<box<P>, box<P>>("BOX(-180 0,-170 1)", "BOX(179 0, 190 1)",  true);
    test_geometry<box<P>, box<P>>("BOX(-180 0,-170 1)", "BOX(181 0, 190 1)",  false); // invalid?
    test_geometry<box<P>, box<P>>("BOX(-180 0,-170 1)", "BOX(180 0, 189 1)",  false); // invalid?

    // Related to https://svn.boost.org/trac/boost/ticket/12412
    test_geometry<box<P>, box<P>>("BOX(-1.346346 51.6, -0.127592 51.7)", "BOX(-2.08882 51.5034, -0.127592 51.9074)", true);
    test_geometry<box<P>, box<P>>("BOX(-2.08882 51.6, -1.346346 51.7)", "BOX(-2.08882 51.5034, -0.127592 51.9074)", true);
    test_geometry<box<P>, box<P>>("BOX(0.127592 51.6, 1.346346 51.7)", "BOX(0.127592 51.5034, 2.08882 51.9074)", true);
    test_geometry<box<P>, box<P>>("BOX(1.346346 51.6, 2.08882 51.7)", "BOX(0.127592 51.5034, 2.08882 51.9074)", true);

    test_geometry<box<P>, box<P>>("BOX(179.08882 1, 180.0 1)", "BOX(179.08882 0, 538.127592 2)", true);
    test_geometry<box<P>, box<P>>("BOX(177.0 1, 178.127592 1)", "BOX(179.08882 0, 538.127592 2)", true);
    test_geometry<box<P>, box<P>>("BOX(179.08882 1, 179.9 1)", "BOX(179.08882 0, 182.127592 2)", true);
    test_geometry<box<P>, box<P>>("BOX(-179.9 1, -177.872408 1)", "BOX(179.08882 0, 182.127592 2)", true);
}

template <typename P>
void test_point_polygon()
{
    std::conditional_t
        <
            std::is_same<typename bg::cs_tag<P>::type, bg::geographic_tag>::value,
            bg::strategy::within::geographic_winding<P>,
            bg::strategy::within::spherical_winding<P>
        > s;

    using poly = bg::model::polygon<P>;

    // MySQL report 08.2017
    test_geometry<P, poly>("POINT(-179 0)",
                           "POLYGON((0 0, 0 2, 2 0, 0 -2, 0 0))",
                           false);
    test_geometry<P, poly>("POINT(-179 0)",
                           "POLYGON((0 0, 0 2, 2 0, 0 -2, 0 0))",
                           false,
                           s);

    test_geometry<P, poly>("POINT(1 0)",
                           "POLYGON((0 0, 0 2, 2 0, 0 -2, 0 0))",
                           true);
    test_geometry<P, poly>("POINT(1 0)",
                           "POLYGON((0 0, 0 2, 2 0, 0 -2, 0 0))",
                           true,
                           s);
}

template <typename P>
void test_cs()
{
    test_point_box<P>();
    test_multi_point_box<P>();
    test_areal_box<P>();
    test_box_box<P>();
    test_point_polygon<P>();
}


int test_main( int , char* [] )
{
    test_cs<bg::model::point<double, 2, bg::cs::spherical_equatorial<bg::degree> > >();
    test_cs<bg::model::point<double, 2, bg::cs::geographic<bg::degree> > >();

    return 0;
}
