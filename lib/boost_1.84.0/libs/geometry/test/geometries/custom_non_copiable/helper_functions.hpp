// Boost.Geometry

// Copyright (c) 2023 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2023 Adam Wulkiewicz, Lodz, Poland.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#ifndef GEOMETRY_TEST_CNC_HELPER_FUNCTIONS
#define GEOMETRY_TEST_CNC_HELPER_FUNCTIONS

#include <boost/geometry.hpp>

#include <sstream>
#include <fstream>

template <typename Ring>
void fill(Ring& ring, std::vector<typename boost::geometry::point_type<Ring>::type> const& v)
{
    ring.custom_clear();
    for(auto const& p : v)
    {
        ring.custom_push_back(p);
    }
}

template <typename Geometry>
auto centroid(Geometry const& geo)
{
    using point_t = typename boost::geometry::point_type<Geometry>::type;
    return boost::geometry::return_centroid<point_t>(geo);
}

template <typename Geometry>
auto envelope(Geometry const& geo)
{
    using point_t = typename boost::geometry::point_type<Geometry>::type;
    using box_t = boost::geometry::model::box<point_t>;
    return boost::geometry::return_envelope<box_t>(geo);
}

template <typename Geometry>
auto hull(Geometry const& geo)
{
    using point_t = typename boost::geometry::point_type<Geometry>::type;
#if defined(BOOST_GEOMETRY_TODO_SUPPORT_CUSTOM_CONVEX_HULL_RESULT)
    // See issue 1136
    using result_t = cnc_polygon<point_t>;
#else
    using result_t = boost::geometry::model::polygon<point_t>;
#endif
    result_t result;
    boost::geometry::convex_hull(geo, result);
    return result;
}

template <typename Geometry>
auto point_on_surface(Geometry const& geo)
{
    using point_t = typename boost::geometry::point_type<Geometry>::type;
    point_t result;
    boost::geometry::point_on_surface(geo, result);
    return result;
}

template <typename ResultGeometry, typename Geometry1, typename Geometry2>
auto areal_intersection(Geometry1 const& a, Geometry2 const& b)
{
    ResultGeometry result;
    boost::geometry::intersection(a, b, result);
    return result;
}

template <typename ResultGeometry, typename Geometry1, typename Geometry2>
auto linear_intersection(Geometry1 const& a, Geometry2 const& b)
{
    ResultGeometry result;
    boost::geometry::intersection(a, b, result);
    return result;
}

template <typename ResultGeometry, typename Geometry>
auto buffer(Geometry const& geo, double distance)
{
    boost::geometry::strategy::buffer::join_round join_strategy(36);
    boost::geometry::strategy::buffer::distance_symmetric<double> distance_strategy(distance);
    boost::geometry::strategy::buffer::end_round end_strategy(36);
    boost::geometry::strategy::buffer::side_straight side_strategy;
    boost::geometry::strategy::buffer::point_circle point_strategy;

    ResultGeometry result;
    boost::geometry::buffer(geo, result,
                            distance_strategy, side_strategy,
                            join_strategy, end_strategy, point_strategy);
    return result;
}

// Currently the "convert" functionality detects if the geometry is the same,
// and if it is, it uses assignment operator. But that does not work for CNC
// Therefore force compilation to another type.
template <typename DestinationType, typename Geometry>
auto convert_to(Geometry const& geo)
{
    DestinationType result;
    boost::geometry::convert(geo, result);
    return result;
}

template <typename Geometry1, typename Geometry2, typename Geometry3>
void create_svg(std::ostream& stream, Geometry1 const& a, Geometry2 const& b, Geometry3 const& intersected)
{
    using point_t = typename boost::geometry::point_type<Geometry1>::type;
    boost::geometry::svg_mapper<point_t> mapper(stream, 1000, 1000);

    mapper.add(a);
    mapper.add(b);

    std::string const prefix
        = boost::geometry::util::is_areal<Geometry1>::value
        ? "opacity:0.4;fill:"
        : "opacity:0.4;stroke-width:10;stroke:";
    mapper.map(a, prefix + "rgb(0,128,0);");
    mapper.map(b, prefix + "rgb(0,0,255);");
    mapper.map(intersected, prefix + "rgb(255,0,0);");
}

void write_svg(std::ostringstream& svg, std::string const& filename)
{
    boost::ignore_unused(svg, filename);
#if defined(TEST_WITH_SVG_FILE)
    std::ofstream tmp("/tmp/" + filename);
    if (tmp.good())
    {
        tmp << svg.str();
    }
#endif

}

#endif // GEOMETRY_TEST_CNC_HELPER_FUNCTIONS
