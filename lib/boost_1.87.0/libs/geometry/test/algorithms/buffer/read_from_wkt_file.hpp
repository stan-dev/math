// Boost.Geometry (aka GGL, Generic Geometry Library)
// Unit Test

// Copyright (c) 2014-2022 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_TEST_READ_FROM_WKT_FILE_HPP
#define BOOST_GEOMETRY_TEST_READ_FROM_WKT_FILE_HPP

#include <boost/geometry/io/wkt/wkt.hpp>

#include <fstream>
#include <string>

// Reads all geometries from WKT file and return it as one big multi polygon WKT
template <typename MultiPolygon>
inline std::string read_from_wkt_file(std::string const& filename)
{
    MultiPolygon mp;
    std::ifstream in(filename.c_str());
    while (in.good())
    {
        std::string line;
        std::getline(in, line);
        if (! line.empty() && line.substr(0, 1) != "#")
        {
            typename boost::range_value<MultiPolygon>::type geometry;
            bg::read_wkt(line, geometry);
            mp.push_back(geometry);
        }
    }
    std::ostringstream out;
    if (! mp.empty())
    {
        out << std::fixed << std::setprecision(20) << bg::wkt(mp);
    }

    BOOST_CHECK(! out.str().empty());

    return out.str();
}

#endif
