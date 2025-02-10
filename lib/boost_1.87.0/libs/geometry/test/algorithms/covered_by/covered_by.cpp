// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2015 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2013-2015 Adam Wulkiewicz, Lodz, Poland.

// This file was modified by Oracle on 2015, 2017, 2022.
// Modifications copyright (c) 2017-2022 Oracle and/or its affiliates.

// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "test_covered_by.hpp"


#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/point_xy.hpp>


template <typename P>
void test_all()
{
    using seg = bg::model::segment<P>;

    test_geometry<P, P>("POINT(0 0)", "POINT(0 0)", true);
    test_geometry<P, P>("POINT(0 0)", "POINT(1 1)", false);

    test_geometry<P, mpt<P>>("POINT(0 0)", "MULTIPOINT(0 0, 1 1)", true);
    test_geometry<mpt<P>, P>("MULTIPOINT(0 0, 1 1)", "POINT(1 1)", false);
    test_geometry<mpt<P>, mpt<P>>("MULTIPOINT(0 0, 1 1)", "MULTIPOINT(1 1, 2 2)", false);


    test_geometry<P, seg>("POINT(1 1)", "LINESTRING(0 0, 2 2)", true);
    test_geometry<P, seg>("POINT(0 0)", "LINESTRING(0 0, 1 1)", true);
    test_geometry<P, seg>("POINT(1 0)", "LINESTRING(0 0, 1 1)", false);

    // linestrings
    test_geometry<P, ls<P>>("POINT(0 0)", "LINESTRING(0 0,1 1,2 2)", true);
    test_geometry<P, ls<P>>("POINT(3 3)", "LINESTRING(0 0,1 1,2 2)", false);
    test_geometry<P, ls<P>>("POINT(1 1)", "LINESTRING(0 0,2 2,3 3)", true);

    // multi_linestrings
    test_geometry<P, mls<P>>("POINT(0 0)", "MULTILINESTRING((0 0,1 1,2 2),(0 0,0 1))", true);
    test_geometry<P, mls<P>>("POINT(0 0)", "MULTILINESTRING((0 0,1 1,2 2),(0 0,0 1),(0 0,1 0))", true);

    // multi_point/segment
    test_geometry<mpt<P>, seg>("MULTIPOINT(0 0, 1 1)", "LINESTRING(0 0, 2 2)", true);

    // multi_point/linestring
    test_geometry<mpt<P>, ls<P>>("MULTIPOINT(0 0, 2 2)", "LINESTRING(0 0, 2 2)", true);
    test_geometry<mpt<P>, ls<P>>("MULTIPOINT(1 1, 3 3)", "LINESTRING(0 0, 2 2)", false);

    // multi_point/multi_linestring
    test_geometry<mpt<P>, mls<P>>("MULTIPOINT(0 0, 1 1)", "MULTILINESTRING((0 0, 2 2),(2 2, 3 3))", true);
    test_geometry<mpt<P>, mls<P>>("MULTIPOINT(0 0, 2 2)", "MULTILINESTRING((0 0, 2 2),(2 2, 3 3))", true);
    test_geometry<mpt<P>, mls<P>>("MULTIPOINT(0 0, 3 3)", "MULTILINESTRING((0 0, 2 2),(2 2, 3 3))", true);
    test_geometry<mpt<P>, mls<P>>("MULTIPOINT(1 1, 4 4)", "MULTILINESTRING((0 0, 2 2),(2 2, 3 3))", false);

    // point/A
    test_geometry<P, ring<P>>("POINT(1 1)", "POLYGON((0 0,0 2,2 2,2 0,0 0))", true);
    test_geometry<P, poly<P>>("POINT(1 1)", "POLYGON((0 0,0 2,2 2,2 0,0 0))", true);
    test_geometry<P, mpoly<P>>("POINT(1 1)", "MULTIPOLYGON(((0 0,0 2,2 2,2 0,0 0)),((2 2,2 3,3 3,3 2,2 2)))", true);

    // on border/corner
    test_geometry<P, poly<P>>("POINT(0 0)", "POLYGON((0 0,0 2,2 2,2 0,0 0))", true);
    test_geometry<P, poly<P>>("POINT(0 1)", "POLYGON((0 0,0 2,2 2,2 0,0 0))", true);

    // aligned to segment/vertex
    test_geometry<P, poly<P>>("POINT(1 1)", "POLYGON((0 0,0 3,3 3,3 1,2 1,2 0,0 0))", true);
    test_geometry<P, poly<P>>("POINT(1 1)", "POLYGON((0 0,0 3,4 3,3 1,2 2,2 0,0 0))", true);

    // same polygon, but point on border
    test_geometry<P, poly<P>>("POINT(3 3)", "POLYGON((0 0,0 3,3 3,3 1,2 1,2 0,0 0))", true);
    test_geometry<P, poly<P>>("POINT(3 3)", "POLYGON((0 0,0 3,4 3,3 1,2 2,2 0,0 0))", true);

    // holes
    test_geometry<P, bg::model::polygon<P> >("POINT(2 2)",
        "POLYGON((0 0,0 4,4 4,4 0,0 0),(1 1,3 1,3 3,1 3,1 1))", false);

    // test multi-with-one-polygon (trivial case)
    test_geometry<P, mpoly<P>>("POINT(1 1)", "MULTIPOLYGON(((0 0,0 2,2 2,2 0,0 0)))", true);
    test_geometry<P, mpoly<P>>("POINT(3 3)", "MULTIPOLYGON(((0 0,0 2,2 2,2 0,0 0)))", false);
    test_geometry<P, mpoly<P>>("POINT(0 1)", "MULTIPOLYGON(((0 0,0 2,2 2,2 0,0 0)))", true);
    test_geometry<P, mpoly<P>>("POINT(4 4)", "MULTIPOLYGON(((0 0,0 2,2 2,2 0,0 0)))", false);

    // test if it is in one of them
    std::string multi("MULTIPOLYGON("
        "((0 0,0 2,2 2,2 0,0 0))"
        "((3 3,3 6,6 6,6 3,3 3))"
        ")");
    test_geometry<P, mpoly<P>>("POINT(4 4)", multi, true);
    test_geometry<P, mpoly<P>>("POINT(1 1)", multi, true);
    test_geometry<P, mpoly<P>>("POINT(0 1)", multi, true);


    // multi_point/A
    test_geometry<mpt<P>, ring<P>>("MULTIPOINT(0 0, 1 1)", "POLYGON((0 0,0 2,2 2,2 0,0 0))", true);
    test_geometry<mpt<P>, poly<P>>("MULTIPOINT(0 0, 2 2)", "POLYGON((0 0,0 2,2 2,2 0,0 0))", true);
    test_geometry<mpt<P>, poly<P>>("MULTIPOINT(1 1, 3 3)", "POLYGON((0 0,0 2,2 2,2 0,0 0))", false);
    test_geometry<mpt<P>, mpoly<P>>("MULTIPOINT(0 0, 1 1)", "MULTIPOLYGON(((0 0,0 2,2 2,2 0,0 0)),((2 2,2 3,3 3,3 2,2 2)))", true);
    test_geometry<mpt<P>, mpoly<P>>("MULTIPOINT(0 0, 2 2)", "MULTIPOLYGON(((0 0,0 2,2 2,2 0,0 0)),((2 2,2 3,3 3,3 2,2 2)))", true);
    test_geometry<mpt<P>, mpoly<P>>("MULTIPOINT(0 0, 3 3)", "MULTIPOLYGON(((0 0,0 2,2 2,2 0,0 0)),((2 2,2 3,3 3,3 2,2 2)))", true);
    test_geometry<mpt<P>, mpoly<P>>("MULTIPOINT(1 1, 4 4)", "MULTIPOLYGON(((0 0,0 2,2 2,2 0,0 0)),((2 2,2 3,3 3,3 2,2 2)))", false);

    test_geometry<P, box<P>>("POINT(1 1)", "BOX(0 0,2 2)", true);
    test_geometry<P, box<P>>("POINT(0 0)", "BOX(0 0,2 2)", true);
    test_geometry<P, box<P>>("POINT(2 2)", "BOX(0 0,2 2)", true);
    test_geometry<P, box<P>>("POINT(0 1)", "BOX(0 0,2 2)", true);
    test_geometry<P, box<P>>("POINT(1 0)", "BOX(0 0,2 2)", true);
    test_geometry<P, box<P>>("POINT(3 3)", "BOX(0 0,2 2)", false);

    test_geometry<mpt<P>, box<P>>("MULTIPOINT(1 1, 2 1)", "BOX(0 0,3 3)", true);
    test_geometry<mpt<P>, box<P>>("MULTIPOINT(0 0, 1 1)", "BOX(0 0,2 2)", true);
    test_geometry<mpt<P>, box<P>>("MULTIPOINT(0 0, 3 4)", "BOX(0 0,2 2)", false);

    test_geometry<ls<P>, box<P>>("LINESTRING(0 0,1 1,1 2)", "BOX(0 0,2 2)", true);
    test_geometry<ls<P>, box<P>>("LINESTRING(0 0,1 1,1 2,1 3)", "BOX(0 0,2 2)", false);
    test_geometry<mls<P>, box<P>>("MULTILINESTRING((0 0,1 1,1 2),(0 1,1 0))", "BOX(0 0,2 2)", true);
    test_geometry<mls<P>, box<P>>("MULTILINESTRING((0 0,1 1,1 2,1 3),(0 1,1 0))", "BOX(0 0,2 2)", false);

    test_geometry<ring<P>, box<P>>("POLYGON((0 0,0 3,3 3,3 0,0 0))", "BOX(0 0,4 4)", true);
    test_geometry<ring<P>, box<P>>("POLYGON((0 0,0 3,3 3,5 0,0 0))", "BOX(0 0,4 4)", false);
    test_geometry<poly<P>, box<P>>("POLYGON((0 0,0 3,3 3,3 0,0 0))", "BOX(0 0,4 4)", true);
    test_geometry<poly<P>, box<P>>("POLYGON((0 0,0 3,3 3,5 0,0 0))", "BOX(0 0,4 4)", false);
    test_geometry<mpoly<P>, box<P>>("MULTIPOLYGON(((0 0,0 3,3 3,3 0,0 0)),((4 4,4 7,7 7,4 7,4 4)))", "BOX(0 0,7 7)", true);
    test_geometry<mpoly<P>, box<P>>("MULTIPOLYGON(((0 0,0 3,3 3,5 0,0 0)),((4 4,4 7,7 7,4 7,4 4)))", "BOX(0 0,4 4)", false);

    test_geometry<box<P>, box<P>>("BOX(1 1,2 2)", "BOX(0 0,3 3)", true);
    test_geometry<box<P>, box<P>>("BOX(0 0,3 3)", "BOX(1 1,2 2)", false);
    test_geometry<box<P>, box<P>>("BOX(0 0,2 2)", "BOX(0 0,3 3)", true);
    test_geometry<box<P>, box<P>>("BOX(1 1,3 3)", "BOX(0 0,3 3)", true);
    test_geometry<box<P>, box<P>>("BOX(1 2,3 3)", "BOX(0 0,3 3)", true);
    test_geometry<box<P>, box<P>>("BOX(1 1,4 3)", "BOX(0 0,3 3)", false);

    test_geometry<box<P>, ring<P>>("BOX(1 1,2 2)", "POLYGON((0 0,0 3,3 3,3 0,0 0))", true);
    test_geometry<box<P>, ring<P>>("BOX(1 1,2 2)", "POLYGON((0 0,0 3,3 1,1 0,0 0))", false);
    test_geometry<box<P>, poly<P>>("BOX(1 1,2 2)", "POLYGON((0 0,0 3,3 3,3 0,0 0))", true);
    test_geometry<box<P>, poly<P>>("BOX(1 1,2 2)", "POLYGON((0 0,0 3,3 1,1 0,0 0))", false);
    test_geometry<box<P>, mpoly<P>>("BOX(1 1,2 2)", "MULTIPOLYGON(((0 0,0 3,3 3,3 0,0 0)),((-1 -1,-3 -4,-7 -7,-4 -3,-1 -1)))", true);
    test_geometry<box<P>, mpoly<P>>("BOX(1 1,2 2)", "MULTIPOLYGON(((0 0,0 3,3 1,1 0,0 0)),((-1 -1,-3 -4,-7 -7,-4 -3,-1 -1)))", false);
}


void test_3d()
{
    using point_type = boost::geometry::model::point<double, 3, boost::geometry::cs::cartesian>;
    box<point_type> box(point_type(0, 0, 0), point_type(4, 4, 4));
    BOOST_CHECK_EQUAL(bg::covered_by(point_type(2, 2, 2), box), true);
    BOOST_CHECK_EQUAL(bg::covered_by(point_type(2, 4, 2), box), true);
    BOOST_CHECK_EQUAL(bg::covered_by(point_type(2, 2, 4), box), true);
    BOOST_CHECK_EQUAL(bg::covered_by(point_type(2, 2, 5), box), false);
}

template <typename P1, typename P2>
void test_mixed_of()
{
    poly<P1> poly1;
    poly<P2> poly2;
    boost::geometry::read_wkt("POLYGON((0 0,0 5,5 5,5 0,0 0))", poly1);
    boost::geometry::read_wkt("POLYGON((0 0,0 5,5 5,5 0,0 0))", poly2);

    box<P1> box1(P1(1, 1), P1(4, 4));
    box<P2> box2(P2(0, 0), P2(5, 5));
    P1 p1(3, 3);
    P2 p2(3, 3);

    BOOST_CHECK_EQUAL(bg::covered_by(p1, poly2), true);
    BOOST_CHECK_EQUAL(bg::covered_by(p2, poly1), true);
    BOOST_CHECK_EQUAL(bg::covered_by(p2, box1), true);
    BOOST_CHECK_EQUAL(bg::covered_by(p1, box2), true);
    BOOST_CHECK_EQUAL(bg::covered_by(box1, box2), true);
    BOOST_CHECK_EQUAL(bg::covered_by(box2, box1), false);

    // TODO: the following does not compile due to incompatible coordinate types
    // (probably needed by overlay)
    //BOOST_CHECK_EQUAL(bg::covered_by(poly1, poly2), true);
    //BOOST_CHECK_EQUAL(bg::covered_by(box1, poly1), true);
    //BOOST_CHECK_EQUAL(bg::covered_by(box2, poly1), true);
}


void test_mixed()
{
    // Mixing point types and coordinate types
    test_mixed_of
        <
            boost::geometry::model::d2::point_xy<double>,
            boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>
        >();
    test_mixed_of
        <
            boost::geometry::model::d2::point_xy<float>,
            boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>
        >();
    test_mixed_of
        <
            boost::geometry::model::d2::point_xy<int>,
            boost::geometry::model::d2::point_xy<double>
        >();
}


int test_main( int , char* [] )
{
    test_all<bg::model::d2::point_xy<int> >();
    test_all<bg::model::d2::point_xy<double> >();

    test_mixed();
    test_3d();

    return 0;
}
