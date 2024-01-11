// Boost.Geometry

// Copyright (c) 2016-2023 Oracle and/or its affiliates.

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
        > ws;

    using poly = bg::model::polygon<P>;

    // MySQL report 08.2017
    test_geometry<P, poly>("POINT(-179 0)",
                           "POLYGON((0 0, 0 2, 2 0, 0 -2, 0 0))",
                           false);
    test_geometry<P, poly>("POINT(-179 0)",
                           "POLYGON((0 0, 0 2, 2 0, 0 -2, 0 0))",
                           false,
                           ws);

    test_geometry<P, poly>("POINT(1 0)",
                           "POLYGON((0 0, 0 2, 2 0, 0 -2, 0 0))",
                           true);
    test_geometry<P, poly>("POINT(1 0)",
                           "POLYGON((0 0, 0 2, 2 0, 0 -2, 0 0))",
                           true,
                           ws);

    using Point = P;
    // Segment going through pole
    {
        bg::model::polygon<Point> poly_n1;
        bg::read_wkt("POLYGON((-90 80,90 80,90 70,-90 70, -90 80))", poly_n1);
        // Points on segment
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-90, 85), poly_n1, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(90, 85), poly_n1, ws), true);
        // Points on pole
        BOOST_CHECK_EQUAL(bg::covered_by(Point(90, 90), poly_n1, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(0, 90), poly_n1, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(45, 90), poly_n1, ws), true);
    }
    // Segment going through pole
    {
        bg::model::polygon<Point> poly_n2;
        bg::read_wkt("POLYGON((-90 80,90 70,0 70,-90 80))", poly_n2);
        // Points on segment
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-90, 85), poly_n2, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(90, 75), poly_n2, ws), true);
        // Points outside but on the same level as segment
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-90, 75), poly_n2, ws), false);
    }
    // Possibly invalid, 2-segment polygon with segment going through pole
    /*{
        bg::model::polygon<Point> poly_n;
        bg::read_wkt("POLYGON((-90 80,90 70,-90 80))", poly_n);
        // Point within
        BOOST_CHECK_EQUAL(bg::within(Point(0, 89), poly_n), true);
        // Points on segment
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-90, 85), poly_n), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(90, 75), poly_n), true);
        // Points outside but on the same level as segment
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-90, 75), poly_n), false);
    }*/
    // Segment endpoints on North pole with arbitrary longitudes
    {
        bg::model::polygon<Point> poly_n4;
        bg::read_wkt("POLYGON((45 90,45 80,-10 80,45 90))", poly_n4);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-10, 85), poly_n4, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(0, 85), poly_n4, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(45, 85), poly_n4, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(30, 85), poly_n4, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(50, 85), poly_n4, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-5, 85), poly_n4, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(30, 70), poly_n4, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(50, 70), poly_n4, ws), false);

        // the same polygon but with two points representing the pole
        bg::model::polygon<Point> poly_n4b;
        bg::read_wkt("POLYGON((45 90,45 80,-10 80,60 90,45 90))", poly_n4b);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-10, 85), poly_n4b, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(0, 85), poly_n4b, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(45, 85), poly_n4b, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(30, 85), poly_n4b, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(50, 85), poly_n4b, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-5, 85), poly_n4b, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(30, 70), poly_n4b, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(50, 70), poly_n4b, ws), false);

        bg::model::polygon<Point> poly_n5;
        bg::read_wkt("POLYGON((0 90,-10 80,45 80,0 90))", poly_n5);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-10, 85), poly_n5, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(0, 85), poly_n5, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-1, 85), poly_n5, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(1, 85), poly_n5, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(45, 85), poly_n5, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(30, 85), poly_n5, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(50, 85), poly_n5, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-5, 85), poly_n5, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(30, 70), poly_n5, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(50, 70), poly_n5, ws), false);

        bg::model::polygon<Point> poly_n_4edges;
        bg::read_wkt("POLYGON((0 90,-10 70,5 60,20 80,0 90))", poly_n_4edges);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(3, 89), poly_n_4edges, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-1, 87), poly_n_4edges, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-10, 86), poly_n_4edges, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(15, 84), poly_n_4edges, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-1, 61), poly_n_4edges, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(15, 81), poly_n_4edges, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(7, 50), poly_n_4edges, ws), false);

        bg::model::polygon<Point> poly_n_5edges;
        bg::read_wkt("POLYGON((0 90,-10 70,5 60,10 85,20 80,0 90))", poly_n_5edges);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(3, 89), poly_n_5edges, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-1, 87), poly_n_5edges, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-10, 86), poly_n_5edges, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(15, 84), poly_n_5edges, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-1, 61), poly_n_5edges, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(15, 81), poly_n_5edges, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(7, 50), poly_n_5edges, ws), false);
    }
    // Segment going through pole
    {
        bg::model::polygon<Point> poly_s1;
        bg::read_wkt("POLYGON((-90 -80,-90 -70,90 -70,90 -80,-90 -80))", poly_s1);
        // Points on segment
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-90, -85), poly_s1, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(90, -85), poly_s1, ws), true);
        // Points on pole
        BOOST_CHECK_EQUAL(bg::covered_by(Point(90, -90), poly_s1, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(0, -90), poly_s1, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(45, -90), poly_s1, ws), true);
    }
    // Segment endpoints on South pole with arbitrary longitudes
    {
        bg::model::polygon<Point> poly_s2;
        bg::read_wkt("POLYGON((45 -90,0 -80,45 -80,45 -90))", poly_s2);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(0, -85), poly_s2, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(45, -85), poly_s2, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(30, -85), poly_s2, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(50, -85), poly_s2, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-5, -85), poly_s2, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(30, -70), poly_s2, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(50, -70), poly_s2, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-5, -70), poly_s2, ws), false);

        bg::model::polygon<Point> poly_s3;
        bg::read_wkt("POLYGON((45 -90,-10 -80,45 -80,45 -90))", poly_s3);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-10, -85), poly_s3, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(0, -85), poly_s3, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-1, -85), poly_s3, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(1, -85), poly_s3, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(45, -85), poly_s3, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(30, -85), poly_s3, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(50, -85), poly_s3, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-5, -85), poly_s3, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(30, -70), poly_s3, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(50, -70), poly_s3, ws), false);

        bg::model::polygon<Point> poly_s5;
        bg::read_wkt("POLYGON((0 -90,-10 -80,45 -80,0 -90))", poly_s5);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-10, -85), poly_s5, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(0, -85), poly_s5, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-1, -85), poly_s5, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(1, -85), poly_s5, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(45, -85), poly_s5, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(30, -85), poly_s5, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(50, -85), poly_s5, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-5, -85), poly_s5, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(30, -70), poly_s5, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(50, -70), poly_s5, ws), false);

        bg::model::polygon<Point> poly_s4;
        bg::read_wkt("POLYGON((0 -89,-10 -80,45 -80,0 -89))", poly_s4);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-10, -85), poly_s4, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(0, -85), poly_s4, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(45, -85), poly_s4, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(30, -85), poly_s4, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(50, -85), poly_s4, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-5, -85), poly_s4, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(30, -71), poly_s4, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(50, -70), poly_s4, ws), false);

        //more complex examples
        bg::model::polygon<Point> poly_s_complex_4edges;
        bg::read_wkt("POLYGON((0 -90,-10 -70,5 -60,20 -80,0 -90))", poly_s_complex_4edges);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(3, -89), poly_s_complex_4edges, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-1, -87), poly_s_complex_4edges, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-10, -86), poly_s_complex_4edges, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(15, -84), poly_s_complex_4edges, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-1, -61), poly_s_complex_4edges, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(15, -81), poly_s_complex_4edges, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(7, -50), poly_s_complex_4edges, ws), false);

        bg::model::polygon<Point> poly_s_complex_5edges;
        bg::read_wkt("POLYGON((0 -90,-10 -70,5 -60,10 -85,20 -80,0 -90))", poly_s_complex_5edges);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(3, -89), poly_s_complex_5edges, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-1, -87), poly_s_complex_5edges, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-10, -86), poly_s_complex_5edges, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(15, -84), poly_s_complex_5edges, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(-1, -61), poly_s_complex_5edges, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(15, -81), poly_s_complex_5edges, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(7, -50), poly_s_complex_5edges, ws), false);
    }
    // Polygon covering nearly half of the globe but no poles
    {
        bg::model::polygon<Point> poly_h1;
        bg::read_wkt("POLYGON((170 0, 170 -80,10 -80,0 -80,0 -20,10 -20,10 20,0 20,0 80,10 80,170 80,170 0))", poly_h1);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, 90), poly_h1, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, 85), poly_h1, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, 50), poly_h1, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, 0), poly_h1, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, -50), poly_h1, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, -85), poly_h1, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, -90), poly_h1, ws), false);
    }
    // Polygon covering more than half of the globe with both holes
    {
        bg::model::polygon<Point> poly_h2;
        bg::read_wkt("POLYGON((180 0, 180 -80,0 -80,10 -80,10 -20,0 -20,0 20,10 20,10 80,0 80,180 80,180 0))", poly_h2);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, 90), poly_h2, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, 85), poly_h2, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, 50), poly_h2, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, 0), poly_h2, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, -50), poly_h2, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, -85), poly_h2, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, -90), poly_h2, ws), true);
    }
    // Polygon covering around half of the globe covering south pole
    {
        bg::model::polygon<Point> poly_h3;
        bg::read_wkt("POLYGON((180 0, 180 -80,0 -80,0 -20,10 -20,10 20,0 20,0 80,10 80,170 80,180 0))", poly_h3);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, 90), poly_h3, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, 85), poly_h3, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, 50), poly_h3, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, 0), poly_h3, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, -50), poly_h3, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, -85), poly_h3, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, -90), poly_h3, ws), true);
    }
    // Polygon covering around half of the globe covering north pole
    {
        bg::model::polygon<Point> poly_h4;
        bg::read_wkt("POLYGON((180 0, 170 -80,10 -80,10 -20,0 -20,0 20,10 20,10 80,0 80,180 80,180 0))", poly_h4);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, 90), poly_h4, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, 85), poly_h4, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, 50), poly_h4, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, 0), poly_h4, ws), true);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, -50), poly_h4, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, -85), poly_h4, ws), false);
        BOOST_CHECK_EQUAL(bg::covered_by(Point(5, -90), poly_h4, ws), false);
    }

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
