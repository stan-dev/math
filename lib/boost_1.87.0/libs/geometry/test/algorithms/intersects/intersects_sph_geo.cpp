// Boost.Geometry

// Copyright (c) 2016 Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "test_intersects.hpp"


#include <boost/geometry/geometries/geometries.hpp>


template <typename P>
void test_point_box()
{
    typedef bg::model::box<P> box_t;

    test_geometry<P, box_t>("POINT(0 0)",    "BOX(0 0, 1 1)", true);
    test_geometry<P, box_t>("POINT(1 1)",    "BOX(0 0, 2 2)", true);

    test_geometry<P, box_t>("POINT(180 1)",  "BOX(170 0, 190 2)", true);
    test_geometry<P, box_t>("POINT(-180 1)", "BOX(170 0, 190 2)", true);
    test_geometry<P, box_t>("POINT(180 1)",  "BOX(170 0, 180 2)", true);
    test_geometry<P, box_t>("POINT(-180 1)", "BOX(170 0, 180 2)", true);
    test_geometry<P, box_t>("POINT(179 1)",  "BOX(170 0, 190 2)", true);
    test_geometry<P, box_t>("POINT(-179 1)", "BOX(170 0, 190 2)", true);
    test_geometry<P, box_t>("POINT(179 1)",  "BOX(170 0, 180 2)", true);
    test_geometry<P, box_t>("POINT(-179 1)", "BOX(170 0, 180 2)", false);
    test_geometry<P, box_t>("POINT(169 1)", "BOX(170 0, 180 2)", false);
}

template <typename P>
void test_box_box()
{
    typedef bg::model::box<P> box_t;

    test_geometry<box_t, box_t>("BOX(0 0, 1 1)", "BOX(0 0, 1 1)", true);

    test_geometry<box_t, box_t>("BOX(-170 0,-160 1)", "BOX(-180 0, 180 1)", true);
    test_geometry<box_t, box_t>("BOX(-170 0,-160 1)", "BOX(170 0, 200 1)",  true);
    test_geometry<box_t, box_t>("BOX(-170 0,-150 1)", "BOX(170 0, 200 1)",  true);
    test_geometry<box_t, box_t>("BOX(201 0,202 1)",   "BOX(170 0, 200 1)",  false); // invalid g1?
    test_geometry<box_t, box_t>("BOX(-159 0,-158 1)", "BOX(170 0, 200 1)",  false);
    test_geometry<box_t, box_t>("BOX(160 0,169 1)",   "BOX(170 0, 200 1)",  false);
    test_geometry<box_t, box_t>("BOX(-159 0,169 1)",  "BOX(170 0, 200 1)",  false);
    test_geometry<box_t, box_t>("BOX(0 0,1 1)",       "BOX(170 0, 370 1)",  true);
    test_geometry<box_t, box_t>("BOX(0 0,10 1)",      "BOX(170 0, 370 1)",  true);
    test_geometry<box_t, box_t>("BOX(-180 0,10 1)",   "BOX(170 0, 370 1)",  true);
    test_geometry<box_t, box_t>("BOX(-180 0,20 1)",   "BOX(170 0, 370 1)",  true);
    test_geometry<box_t, box_t>("BOX(10 0,20 1)",     "BOX(170 0, 370 1)",  true);
    test_geometry<box_t, box_t>("BOX(160 0,180 1)",   "BOX(170 0, 370 1)",  true);
    test_geometry<box_t, box_t>("BOX(160 0,165 1)",   "BOX(170 0, 370 1)",  false);
    test_geometry<box_t, box_t>("BOX(15 0,20 1)",     "BOX(170 0, 370 1)",  false);
    test_geometry<box_t, box_t>("BOX(375 0,380 1)",   "BOX(170 0, 370 1)",  false); // invalid g1?

    test_geometry<box_t, box_t>("BOX(-180 0,-170 1)", "BOX(180 0, 190 1)",  true); // invalid?
    test_geometry<box_t, box_t>("BOX(-180 0,-170 1)", "BOX(180 0, 191 1)",  true); // invalid?
    test_geometry<box_t, box_t>("BOX(-180 0,-170 1)", "BOX(179 0, 190 1)",  true);
    test_geometry<box_t, box_t>("BOX(-180 0,-170 1)", "BOX(181 0, 190 1)",  true); // invalid?
    test_geometry<box_t, box_t>("BOX(-180 0,-170 1)", "BOX(180 0, 189 1)",  true); // invalid?
}


template <typename P>
void test_cs()
{
    test_point_box<P>();
    test_box_box<P>();
}


int test_main( int , char* [] )
{
    test_cs<bg::model::point<double, 2, bg::cs::spherical_equatorial<bg::degree> > >();
    test_cs<bg::model::point<double, 2, bg::cs::geographic<bg::degree> > >();

    return 0;
}

//https://github.com/boostorg/geometry/issues/480
namespace spherical {
  using point = boost::geometry::model::d2::point_xy<double, boost::geometry::cs::spherical_equatorial<boost::geometry::degree> >;
  using linestring = boost::geometry::model::linestring<point>;
  using polygon = boost::geometry::model::polygon<point>;
  using multipolygon = boost::geometry::model::multi_polygon<polygon>;
  using mpoint = boost::geometry::model::multi_point<point>;
}

BOOST_AUTO_TEST_CASE(usaIntersectBigTile90) {
  auto const big_tile_wkt = "POLYGON((-182.25 -2.25,-182.25 85,-90 85,4.5 85,4.5 -2.25,-90 -2.25,-182.25 -2.25))";
  auto const big_tile_wkt_on90 = "POLYGON((-182.25 -2.25,-182.25 90,-90 90,4.5 90,4.5 -2.25,-90 -2.25,-182.25 -2.25))";
  auto const usa_wkt = "POLYGON((-125.0840938165823 48.18451089135078,-123.2484442951568 48.28402780487058,-123.3222396381058 49.00207166780857,-106.2294164816153 48.99936106028161,-94.95738881549684 49.37019439405368,-91.56749988130491 48.04377777876205,-88.36986182651864 48.30606296052222,-82.51861102256366 45.33863888434129,-82.67972207213995 41.67655556472526,-74.97255548445517 44.98347221572623,-71.50101930291031 45.01344579588383,-69.22446981284281 47.4598397141629,-67.79019834014535 47.06722478688557,-67.04184250949656 44.36199518537477,-70.23005912341451 43.2124078816444,-69.75267183371977 41.13341339615248,-73.71086681700677 40.38246941502378,-75.73925419455345 36.92336928786909,-75.32570262625269 35.10289000177789,-81.00038138644489 31.36986317745486,-80.21558018606714 24.87617100882401,-82.31642199133823 24.40399007010678,-81.37277484337314 25.35110770766313,-82.98991779576328 28.87424738341057,-86.63551567900944 30.19587268651689,-88.60142307258893 30.01998991835538,-89.40197910608261 28.705405956028,-94.4542152059084 29.30489831467678,-96.8117767445914 27.80056771982362,-97.42440989188675 25.84041996728249,-99.0845302237628 26.39797732094523,-101.4011198718317 29.77120000013631,-103.2862798604784 28.97691995351909,-106.4506798904261 31.76422996252971,-111.074824919583 31.33223902311842,-114.7199604399525 32.71865525370477,-117.4636899326359 32.58936795833381,-118.1537433047147 33.49758670948682,-118.5200689113933 32.61376470726387,-119.9333407371885 33.34100477748434,-118.6596572250774 33.83057945961541,-120.7240783642898 33.9956452157932,-124.7155563780111 40.37995067155807,-125.0840938165823 48.18451089135078))";

  spherical::polygon big_tile, big_tile_on90, usa;
  boost::geometry::read_wkt(big_tile_wkt, big_tile);
  boost::geometry::read_wkt(big_tile_wkt_on90, big_tile_on90);
  boost::geometry::read_wkt(usa_wkt, usa);

  spherical::multipolygon out1, out2;
  boost::geometry::intersection(big_tile_on90, usa, out1);
  boost::geometry::intersection(big_tile, usa, out2);
  BOOST_CHECK(out1.size() == 1);
  BOOST_CHECK(boost::geometry::area(out1) == boost::geometry::area(out2));
}
