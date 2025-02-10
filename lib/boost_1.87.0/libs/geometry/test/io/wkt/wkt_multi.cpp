// Boost.Geometry (aka GGL, Generic Geometry Library)
// Unit Test

// Copyright (c) 2007-2022 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2015 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2015 Mateusz Loskot, London, UK.

// This file was modified by Oracle on 2014, 2015.
// Modifications copyright (c) 2014-2015 Oracle and/or its affiliates.

// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//#include <iostream>
#include <sstream>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/concept/requires.hpp>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/included/test_exec_monitor.hpp>

#include <boost/geometry/geometries/geometries.hpp>

#include <boost/geometry/io/wkt/wkt.hpp>

template <typename T>
void test_all();


// Include the single test
#define GEOMETRY_TEST_MULTI
#include "io/wkt/wkt.cpp"

template <typename T>
void test_order_closure()
{
    using namespace boost::geometry;
    typedef bg::model::point<T, 2, bg::cs::cartesian> Pt;
    typedef bg::model::polygon<Pt, true, true> PCWC;
    typedef bg::model::polygon<Pt, true, false> PCWO;
    typedef bg::model::polygon<Pt, false, true> PCCWC;
    typedef bg::model::polygon<Pt, false, false> PCCWO;
    typedef bg::model::multi_polygon<PCWC> MPCWC;
    typedef bg::model::multi_polygon<PCWO> MPCWO;
    typedef bg::model::multi_polygon<PCCWC> MPCCWC;
    typedef bg::model::multi_polygon<PCCWO> MPCCWO;

    std::string wkt_cwc = "MULTIPOLYGON(((0 0,0 2,2 2,2 0,0 0)),((0 0,0 -3,-3 -3,-3 0,0 0),(-1 -1,-2 -1,-2 -2,-1 -2,-1 -1)))";
    std::string wkt_cwo = "MULTIPOLYGON(((0 0,0 2,2 2,2 0)),((0 0,0 -3,-3 -3,-3 0),(-1 -1,-2 -1,-2 -2,-1 -2)))";
    std::string wkt_ccwc = "MULTIPOLYGON(((0 0,2 0,2 2,0 2,0 0)),((0 0,-3 0,-3 -3,0 -3,0 0),(-1 -1,-1 -2,-2 -2,-2 -1,-1 -1)))";
    std::string wkt_ccwo = "MULTIPOLYGON(((0 0,2 0,2 2,0 2)),((0 0,-3 0,-3 -3,0 -3),(-1 -1,-1 -2,-2 -2,-2 -1)))";

    test_wkt<MPCWC>(wkt_cwc, wkt_cwc, 15, 0, 12, 24);
    test_wkt<MPCWO>(wkt_cwc, wkt_cwc, 12, 0, 12, 24);
    test_wkt<MPCWO>(wkt_cwo, wkt_cwc, 12, 0, 12, 24);
    test_wkt<MPCCWC>(wkt_ccwc, wkt_ccwc, 15, 0, 12, 24);
    test_wkt<MPCCWO>(wkt_ccwc, wkt_ccwc, 12, 0, 12, 24);
    test_wkt<MPCCWO>(wkt_ccwo, wkt_ccwc, 12, 0, 12, 24);
}

template <typename T>
void test_all()
{
    using namespace boost::geometry;
    typedef bg::model::point<T, 2, bg::cs::cartesian> P;

    test_wkt<bg::model::multi_point<P> >("multipoint((1 2),(3 4))", 2);
    test_wkt<bg::model::multi_linestring<bg::model::linestring<P> > >("multilinestring((1 1,2 2,3 3),(4 4,5 5,6 6))", 6, 4 * sqrt(2.0));
    test_wkt<bg::model::multi_polygon<bg::model::polygon<P> > >("multipolygon(((0 0,0 2,2 2,2 0,0 0),(1 1,1 2,2 2,2 1,1 1)),((0 0,0 4,4 4,4 0,0 0)))", 15, 0, 21, 28);

    // Support tabs, and new lines.
    test_relaxed_wkt<bg::model::multi_point<P> >("multipoint((1\t2),\n(3\r4))", "multipoint((1 2),(3 4))");

    // Verify empty multi geometries
    test_relaxed_wkt<bg::model::multi_point<P> >("multipoint()", "multipoint()");
    test_relaxed_wkt<bg::model::multi_linestring<bg::model::linestring<P> >>("multilinestring()", "multilinestring()");
    test_relaxed_wkt<bg::model::multi_polygon<bg::model::polygon<P> > >("multipolygon()", "multipolygon()");

    test_relaxed_wkt<bg::model::multi_point<P> >("multipoint empty", "multipoint()");
    test_relaxed_wkt<bg::model::multi_linestring<bg::model::linestring<P> >>("multilinestring empty", "multilinestring()");
    test_relaxed_wkt<bg::model::multi_polygon<bg::model::polygon<P> > >("multipolygon empty", "multipolygon()");

    // Support for the official alternative syntax for multipoint
    // (provided by Aleksey Tulinov):
    test_relaxed_wkt<bg::model::multi_point<P> >("multipoint(1 2,3 4)", "multipoint((1 2),(3 4))");

    test_wrong_wkt<bg::model::multi_polygon<bg::model::polygon<P> > >(
        "MULTIPOLYGON(((0 0,0 2,2 2,2 0,0 0),(1 1,1 2,2 2,2 1,1 1)),(0 0,0 4,4 4,4 0,0 0)))",
        "expected '('");

    test_wrong_wkt<bg::model::multi_linestring<bg::model::linestring<P> > >(
        "MULTILINESTRING ((10 10, 20 20, 10 40), (40 40, 30 30, 40 20, 30 10)), (0 0, 1 1)",
        "too many tokens at ','");

    test_wrong_wkt<bg::model::multi_point<P> >(
        "MULTIPOINT((8 9), 10 11)",
        "expected '(' at '10'");
    test_wrong_wkt<bg::model::multi_point<P> >(
        "MULTIPOINT(12 13, (14 15))",
        "bad lexical cast: source type value could not be interpreted as target at '(' in 'multipoint(12 13, (14 15))'");
    test_wrong_wkt<bg::model::multi_point<P> >(
        "MULTIPOINT((16 17), (18 19)",
        "expected ')' in 'multipoint((16 17), (18 19)'");
    test_wrong_wkt<bg::model::multi_point<P> >(
        "MULTIPOINT(16 17), (18 19)",
        "too many tokens at ',' in 'multipoint(16 17), (18 19)'");

    test_order_closure<T>();
}
