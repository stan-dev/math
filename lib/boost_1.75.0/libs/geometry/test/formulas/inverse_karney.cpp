// Boost.Geometry
// Unit Test

// Copyright (c) 2016-2019 Oracle and/or its affiliates.

// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Copyright (c) 2018 Adeel Ahmad, Islamabad, Pakistan.

// Contributed and/or modified by Adeel Ahmad, as part of Google Summer of Code 2018 program

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#include <sstream>

#include "test_formula.hpp"
#include "inverse_cases.hpp"
#include "inverse_cases_antipodal.hpp"
#include "inverse_cases_small_angles.hpp"

#include <boost/geometry/formulas/karney_inverse.hpp>

#include <boost/geometry/srs/spheroid.hpp>

template <typename Result>
void check_inverse(std::string const& name,
                   Result const& results,
                   bg::formula::result_inverse<double> const& result,
                   expected_result const& expected,
                   expected_result const& reference,
                   double reference_error)
{
    std::stringstream ss;
    ss << "(" << results.p1.lon << " " << results.p1.lat << ")->(" << results.p2.lon << " " << results.p2.lat << ")";

    check_one(name + "_d  " + ss.str(),
              result.distance, expected.distance, reference.distance, reference_error);
    check_one(name + "_a  " + ss.str(),
              result.azimuth, expected.azimuth, reference.azimuth, reference_error, true);
    check_one(name + "_ra " + ss.str(),
              result.reverse_azimuth, expected.reverse_azimuth, reference.reverse_azimuth, reference_error, true);
    check_one(name + "_rl " + ss.str(),
              result.reduced_length, expected.reduced_length, reference.reduced_length, reference_error);
    check_one(name + "_gs " + ss.str(),
              result.geodesic_scale, expected.geodesic_scale, reference.geodesic_scale, reference_error);
}

void test_all(expected_results const& results)
{
    double lon1d = results.p1.lon;
    double lat1d = results.p1.lat;
    double lon2d = results.p2.lon;
    double lat2d = results.p2.lat;

    // WGS84
    bg::srs::spheroid<double> spheroid(6378137.0, 6356752.3142451793);

    bg::formula::result_inverse<double> result_k;

    typedef bg::formula::karney_inverse<double, true, true, true, true, true, 8> ka_t;
    result_k = ka_t::apply(lon1d, lat1d, lon2d, lat2d, spheroid);
    check_inverse("karney", results, result_k, results.vincenty, results.reference, 0.0000001);
}

template <typename ExpectedResults>
void test_karney(ExpectedResults const& results)
{
    double lon1d = results.p1.lon;
    double lat1d = results.p1.lat;
    double lon2d = results.p2.lon;
    double lat2d = results.p2.lat;

    // WGS84
    bg::srs::spheroid<double> spheroid(6378137.0, 6356752.3142451793);

    bg::formula::result_inverse<double> result;

    typedef bg::formula::karney_inverse<double, true, true, true, true, true, 8> ka_t;
    result = ka_t::apply(lon1d, lat1d, lon2d, lat2d, spheroid);
    check_inverse("karney", results, result, results.karney, results.karney, 0.0000001);
}

int test_main(int, char*[])
{
    for (size_t i = 0; i < expected_size; ++i)
    {
        test_all(expected[i]);
    }

    for (size_t i = 0; i < expected_size_antipodal; ++i)
    {
        test_karney(expected_antipodal[i]);
    }

    for (size_t i = 0; i < expected_size_small_angles; ++i)
    {
        test_karney(expected_small_angles[i]);
    }

    return 0;
}
