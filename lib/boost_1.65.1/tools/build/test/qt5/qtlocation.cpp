// (c) Copyright Juergen Hunold 2012
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE QtPositioning

#include <QGeoAddress>
#include <QGeoLocation>

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(defines) {
  BOOST_CHECK_EQUAL(BOOST_IS_DEFINED(QT_CORE_LIB), true);
  BOOST_CHECK_EQUAL(BOOST_IS_DEFINED(QT_POSITIONING_LIB), true);
  BOOST_CHECK_EQUAL(BOOST_IS_DEFINED(QT_NETWORK_LIB), true);
  BOOST_CHECK_EQUAL(BOOST_IS_DEFINED(QT_LOCATION_LIB), true);
}

BOOST_TEST_DONT_PRINT_LOG_VALUE(QGeoAddress)

BOOST_AUTO_TEST_CASE(geo_location) {
  QGeoLocation geolocation;

  QGeoAddress address;

  BOOST_CHECK_EQUAL(geolocation.address(), address);
}
