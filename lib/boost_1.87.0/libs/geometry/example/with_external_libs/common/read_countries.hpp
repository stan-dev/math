// Boost.Geometry
//
// Copyright (c) 2007-2024 Barend Gehrels, Amsterdam, the Netherlands.
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef READ_COUNTRIES_HPP
#define READ_COUNTRIES_HPP

#include <fstream>
#include <boost/geometry/geometry.hpp>

// ----------------------------------------------------------------------------
// Read an ASCII file containing WKT's of either POLYGON or MULTIPOLYGON
// ----------------------------------------------------------------------------
template <typename Geometry>
std::vector<Geometry> read_countries(std::string const& filename)
{
    std::vector<Geometry> geometries;
    std::ifstream cpp_file(filename.c_str());
    if (!cpp_file.is_open())
    {
        return geometries;
    }
    while (! cpp_file.eof() )
    {
        std::string line;
        std::getline(cpp_file, line);
        if (line.empty())
        {
            continue;
        }
        Geometry geometry;
        if (line.substr(0, 4) == "POLY")
        {
            using polygon_t = std::decay_t<decltype(*geometry.begin())>;
            polygon_t polygon;
            boost::geometry::read_wkt(line, polygon);
            geometry.push_back(polygon);
        }
        else
        {
            boost::geometry::read_wkt(line, geometry);
        }

        geometries.push_back(geometry);
    }
    return geometries;
}

// Returns the envelope of a collection of geometries
template <typename Box, typename Countries>
Box calculate_envelope(Countries const& countries)
{
    Box box;
    boost::geometry::assign_inverse(box);
    
    for (auto const& country : countries)
    {
        boost::geometry::expand(box, boost::geometry::return_envelope<Box>(country));
    }
    return box;
}

#endif
