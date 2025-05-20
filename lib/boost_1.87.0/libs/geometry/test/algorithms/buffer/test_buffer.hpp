// Boost.Geometry (aka GGL, Generic Geometry Library)
// Unit Test Helper

// Copyright (c) 2010-2019 Barend Gehrels, Amsterdam, the Netherlands.

// This file was modified by Oracle on 2016-2024.
// Modifications copyright (c) 2016-2024, Oracle and/or its affiliates.
// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_GEOMETRY_TEST_BUFFER_HPP
#define BOOST_GEOMETRY_TEST_BUFFER_HPP

#include <iostream>
#include <fstream>
#include <iomanip>

#include <geometry_test_common.hpp>
#include <expectation_limits.hpp>

#include <boost/geometry/algorithms/envelope.hpp>
#include <boost/geometry/algorithms/area.hpp>
#include <boost/geometry/algorithms/buffer.hpp>
#include <boost/geometry/algorithms/correct.hpp>
#include <boost/geometry/algorithms/disjoint.hpp>
#include <boost/geometry/algorithms/intersects.hpp>
#include <boost/geometry/algorithms/is_empty.hpp>
#include <boost/geometry/algorithms/is_valid.hpp>
#include <boost/geometry/algorithms/num_interior_rings.hpp>
#include <boost/geometry/algorithms/union.hpp>

#include <boost/geometry/algorithms/detail/overlay/debug_turn_info.hpp>
#include <boost/geometry/algorithms/detail/overlay/self_turn_points.hpp>

#include <boost/geometry/geometries/geometries.hpp>

#include <boost/geometry/strategies/strategies.hpp>

#include <boost/geometry/strategies/buffer.hpp>

#include <boost/geometry/io/wkt/wkt.hpp>

#include <boost/geometry/util/condition.hpp>

const double same_distance = -999;

#if defined(TEST_WITH_SVG)
#  include "test_buffer_svg.hpp"
#  include "test_buffer_svg_per_turn.hpp"
#endif

#if defined(TEST_WITH_CSV)
#  include "test_buffer_csv.hpp"
#endif

//-----------------------------------------------------------------------------
template <typename JoinStrategy>
struct JoinTestProperties
{
    static std::string name() { return "join_unknown"; }
    static bool is_miter() { return false; }
};

template<> struct JoinTestProperties<boost::geometry::strategy::buffer::join_round>
{
    static std::string name() { return "round"; }
    static bool is_miter() { return false; }
};

template<typename F, typename S, typename CT>
struct JoinTestProperties<boost::geometry::strategy::buffer::geographic_join_round<F, S, CT> >
{
    static std::string name() { return "geo_round"; }
    static bool is_miter() { return false; }
};

template<> struct JoinTestProperties<boost::geometry::strategy::buffer::join_miter>
{
    static std::string name() { return "miter"; }
    static bool is_miter() { return true; }
};

template<typename F, typename S, typename CT>
struct JoinTestProperties<boost::geometry::strategy::buffer::geographic_join_miter<F, S, CT> >
{
    static std::string name() { return "geo_miter"; }
    static bool is_miter() { return true; }
};
template<> struct JoinTestProperties<boost::geometry::strategy::buffer::join_round_by_divide>
{
    static std::string name() { return "divide"; }
    static bool is_miter() { return false; }
};


//-----------------------------------------------------------------------------
template <typename EndStrategy>
struct EndTestProperties { };

template<> struct EndTestProperties<boost::geometry::strategy::buffer::end_round>
{
    static std::string name() { return "round"; }
    static bool is_round() { return true; }
};

template<typename F, typename S, typename CT>
struct EndTestProperties<boost::geometry::strategy::buffer::geographic_end_round<F, S, CT>>
{
    static std::string name() { return "geo_round"; }
    static bool is_round() { return true; }
};

template<> struct EndTestProperties<boost::geometry::strategy::buffer::end_flat>
{
    static std::string name() { return "flat"; }
    static bool is_round() { return false; }
};

struct ut_settings : public ut_base_settings
{
    explicit ut_settings(double tol = 0.01, bool val = true, int points = 88)
        : ut_base_settings(val)
        , tolerance(tol)
        , points_per_circle(points)
    {}

    static inline ut_settings ignore_validity()
    {
        ut_settings result;
        result.set_test_validity(false);
        return result;
    }

    static inline ut_settings assertions_only()
    {
        ut_settings result;
        result.test_area = false;
        result.set_test_validity(false);
        return result;
    }

    static inline double ignore_area() { return 9999.9; }

    double tolerance;
    bool test_area = true;
    bool use_ln_area = false;

    // Number of points in a circle. Not used for geo tests.
    int points_per_circle;

    double multiplier_min_area = 0.95;
    double multiplier_max_area = 1.05;
    double fraction_buffered_points_too_close = 0.10;
};

template
<
    typename GeometryOut,
    typename JoinStrategy,
    typename EndStrategy,
    typename DistanceStrategy,
    typename SideStrategy,
    typename PointStrategy,
    typename Strategy,
    typename Geometry
>
void test_buffer(std::string const& caseid,
            bg::model::multi_polygon<GeometryOut>& buffered,
            Geometry const& geometry,
            JoinStrategy const& join_strategy,
            EndStrategy const& end_strategy,
            DistanceStrategy const& distance_strategy,
            SideStrategy const& side_strategy,
            PointStrategy const& point_strategy,
            Strategy const& strategy,
            int expected_count,
            int expected_holes_count,
            expectation_limits const& expected_area,
            ut_settings const& settings)
{
    namespace bg = boost::geometry;

    typedef typename bg::coordinate_type<Geometry>::type coordinate_type;
    typedef typename bg::point_type<Geometry>::type point_type;

    typedef typename bg::tag<Geometry>::type tag;
    // TODO use something different here:
    std::string type = std::is_same<tag, bg::polygon_tag>::value ? "poly"
        : std::is_same<tag, bg::linestring_tag>::value ? "line"
        : std::is_same<tag, bg::point_tag>::value ? "point"
        : std::is_same<tag, bg::multi_polygon_tag>::value ? "multipoly"
        : std::is_same<tag, bg::multi_linestring_tag>::value ? "multiline"
        : std::is_same<tag, bg::multi_point_tag>::value ? "multipoint"
        : ""
        ;

    bg::model::box<point_type> envelope;
    if (bg::is_empty(geometry))
    {
        bg::assign_values(envelope, 0, 0, 1,  1);
    }
    else
    {
        bg::envelope(geometry, envelope, strategy);
    }

    std::string join_name = JoinTestProperties<JoinStrategy>::name();
    std::string end_name = EndTestProperties<EndStrategy>::name();

    if ( BOOST_GEOMETRY_CONDITION((
            std::is_same<tag, bg::point_tag>::value
         || std::is_same<tag, bg::multi_point_tag>::value )) )
    {
        join_name.clear();
    }

    std::ostringstream complete;
    complete
        << type << "_"
        << caseid << "_"
        << string_from_type<coordinate_type>::name()
        << "_" << join_name
        << (end_name.empty() ? "" : "_") << end_name
        << (distance_strategy.negative() ? "_deflate" : "")
        << (bg::point_order<GeometryOut>::value == bg::counterclockwise ? "_ccw" : "")
         // << "_" << point_buffer_count
        ;

    //std::cout << complete.str() << std::endl;

#if defined(TEST_WITH_CSV)
    detail::buffer_visitor_csv visitor("/tmp/csv/" + caseid + "_");

#elif defined(TEST_WITH_SVG_PER_TURN)
    save_turns_visitor<point_type> visitor;
#elif defined(TEST_WITH_SVG)

    buffer_svg_mapper<point_type> buffer_mapper(complete.str());

    std::ostringstream filename;
    filename << "buffer_" << complete.str() << ".svg";
    std::ofstream svg(filename.str().c_str());
    typedef bg::svg_mapper<point_type> mapper_type;
    mapper_type mapper(svg, 1000, 800);

    svg_visitor<mapper_type, bg::model::box<point_type> > visitor(mapper);

    // Set the SVG boundingbox, with a margin. The margin is necessary because
    // drawing is already started before the buffer is finished. It is not
    // possible to "add" the buffer (unless we buffer twice).
    double margin = distance_strategy.negative()
            ? 1.0
            : 1.1 * distance_strategy.max_distance(join_strategy, end_strategy);

    if (std::is_same<typename bg::coordinate_system<point_type>::type, bg::cs::geographic<bg::degree> >::value)
    {
        // Divide to avoid a too zoomed out SVG.
        // TODO: this can go if bg::buffer for box accepts geographic boxes.
        margin *= 1.25e-5;
    }

    buffer_mapper.prepare(mapper, visitor, envelope, margin);
#else
    bg::detail::buffer::visit_pieces_default_policy visitor;
#endif

    typedef typename bg::point_type<Geometry>::type point_type;

    bg::buffer(envelope, envelope, distance_strategy.max_distance(join_strategy, end_strategy));

    buffered.clear();
    bg::detail::buffer::buffer_inserter<GeometryOut>(geometry,
                        std::back_inserter(buffered),
                        distance_strategy,
                        side_strategy,
                        join_strategy,
                        end_strategy,
                        point_strategy,
                        strategy,
                        visitor);

#if defined(TEST_WITH_CSV)
#elif defined(TEST_WITH_SVG)
    buffer_mapper.map_input_output(mapper, geometry, buffered, distance_strategy.negative());
#endif

    //Uncomment to create simple CSV to compare/use in tests - adapt precision if necessary
    //std::cout << complete.str() << "," << std::fixed << std::setprecision(0) << area << std::endl;
    //return;

    if (bg::is_empty(buffered) && expected_area.is_zero())
    {
        // As expected - don't get rescale policy for output (will be invalid)
        return;
    }

    if (settings.test_area)
    {
        BOOST_CHECK_MESSAGE
            (
                ! bg::is_empty(buffered),
                complete.str() << " output is empty (unexpected)."
            );
    }

    bg::model::box<point_type> envelope_output;
    bg::assign_values(envelope_output, 0, 0, 1,  1);
    bg::envelope(buffered, envelope_output, strategy);

    //    std::cout << caseid << std::endl;
    //    std::cout << "INPUT: " << bg::wkt(geometry) << std::endl;
    //    std::cout << "OUTPUT: " << area << std::endl;
    //    std::cout << "OUTPUT env: " << bg::wkt(envelope_output) << std::endl;
    //    std::cout << bg::wkt(buffered) << std::endl;

    if (expected_count >= 0)
    {
        BOOST_CHECK_MESSAGE
            (
                int(buffered.size()) == expected_count,
                "#outputs not as expected."
                << " Expected: " << expected_count
                << " Detected: " << buffered.size()
            );
    }

    if (expected_holes_count >= 0)
    {
        std::size_t nholes = bg::num_interior_rings(buffered);
        BOOST_CHECK_MESSAGE
        (
            int(nholes) == expected_holes_count,
            complete.str() << " #holes not as expected."
            << " Expected: " << expected_holes_count
            << " Detected: " << nholes
        );
    }

    if (settings.test_area)
    {
        auto const area = bg::area(buffered, strategy);
        BOOST_CHECK_MESSAGE(expected_area.contains(area, settings.tolerance, settings.use_ln_area),
              "difference: " << caseid << std::setprecision(20)
              << " #area expected: " << expected_area
              << " detected: " << area
              << " type: " << (type_for_assert_message<Geometry, GeometryOut>())
              );
    }

    if (settings.test_validity() && ! bg::is_valid(buffered))
    {
        BOOST_CHECK_MESSAGE(bg::is_valid(buffered), complete.str() <<  " is not valid");
    }

#if defined(TEST_WITH_CSV)
    visitor.write_input_output(geometry, buffered);
#elif defined(TEST_WITH_SVG_PER_TURN)
    {
        // Create a per turn visitor to map per turn, and buffer again with it
        per_turn_visitor<point_type> ptv(complete.str(), visitor.get_points());
        bg::detail::buffer::buffer_inserter<GeometryOut>(geometry,
                            std::back_inserter(buffered),
                            distance_strategy,
                            side_strategy,
                            join_strategy,
                            end_strategy,
                            point_strategy,
                            ptv);
        ptv.map_input_output(geometry, buffered, distance_strategy.negative());
        // self_ips NYI here
    }
#elif defined(TEST_WITH_SVG)
    buffer_mapper.map_self_ips(mapper, buffered, strategy);
#endif

}

template
<
    typename GeometryOut,
    typename JoinStrategy,
    typename EndStrategy,
    typename DistanceStrategy,
    typename SideStrategy,
    typename PointStrategy,
    typename Geometry
>
void test_buffer(std::string const& caseid, bg::model::multi_polygon<GeometryOut>& buffered, Geometry const& geometry,
            JoinStrategy const& join_strategy,
            EndStrategy const& end_strategy,
            DistanceStrategy const& distance_strategy,
            SideStrategy const& side_strategy,
            PointStrategy const& point_strategy,
            expectation_limits const& expected_area,
            ut_settings const& settings = ut_settings())
{
    typename bg::strategies::buffer::services::default_strategy
        <
            Geometry
        >::type strategies;

    test_buffer<GeometryOut>(caseid, buffered, geometry,
        join_strategy, end_strategy, distance_strategy, side_strategy, point_strategy,
        strategies,
        -1, -1, expected_area, settings);
}

#ifdef BOOST_GEOMETRY_CHECK_WITH_POSTGIS
static int counter = 0;
#endif

template
<
    typename Geometry,
    typename GeometryOut,
    typename JoinStrategy,
    typename EndStrategy
>
void test_one(std::string const& caseid, std::string const& wkt,
        JoinStrategy const& join_strategy, EndStrategy const& end_strategy,
        int expected_count, int expected_holes_count, expectation_limits const& expected_area,
        double distance_left, ut_settings const& settings = ut_settings(),
        double distance_right = same_distance)
{
    namespace bg = boost::geometry;
    Geometry g;
    bg::read_wkt(wkt, g);
    bg::correct(g);

#ifdef BOOST_GEOMETRY_CHECK_WITH_POSTGIS
    std::cout
        << (counter > 0 ? "union " : "")
        << "select " << counter++
        << ", '" << caseid << "' as caseid"
        << ", ST_Area(ST_Buffer(ST_GeomFromText('" << wkt << "'), "
        << distance_left
        << ", 'endcap=" << end_name << " join=" << join_name << "'))"
        << ", "  << expected_area
        << std::endl;
#endif


    bg::strategy::buffer::side_straight side_strategy;
    bg::strategy::buffer::point_circle circle_strategy(settings.points_per_circle);

    bg::strategy::buffer::distance_asymmetric
    <
        typename bg::coordinate_type<Geometry>::type
    > distance_strategy(distance_left,
                        bg::math::equals(distance_right, same_distance)
                        ? distance_left : distance_right);

    typename bg::strategies::buffer::services::default_strategy
        <
            Geometry
        >::type strategies;

    bg::model::multi_polygon<GeometryOut> buffered;
    test_buffer<GeometryOut>
            (caseid, buffered, g,
            join_strategy, end_strategy,
            distance_strategy, side_strategy, circle_strategy,
            strategies,
            expected_count, expected_holes_count, expected_area,
            settings);

#if !defined(BOOST_GEOMETRY_COMPILER_MODE_DEBUG) \
    && !defined(BOOST_GEOMETRY_TEST_ONLY_ONE_ORDER) \
    && defined(BOOST_GEOMETRY_COMPILER_MODE_RELEASE)

    // Also test symmetric distance strategy if right-distance is not specified
    // (only in release mode, not if "one order" if speficied)
    if (bg::math::equals(distance_right, same_distance))
    {
        bg::strategy::buffer::distance_symmetric
        <
            typename bg::coordinate_type<Geometry>::type
        > sym_distance_strategy(distance_left);

        test_buffer<GeometryOut>
                (caseid + "_sym", buffered, g,
                join_strategy, end_strategy,
                sym_distance_strategy, side_strategy, circle_strategy,
                strategies,
                expected_count, expected_holes_count, expected_area,
                settings);

    }
#endif
}

template
<
    typename Geometry,
    typename GeometryOut,
    typename JoinStrategy,
    typename EndStrategy
>
void test_one(std::string const& caseid, std::string const& wkt,
        JoinStrategy const& join_strategy, EndStrategy const& end_strategy,
        expectation_limits const& expected_area,
        double distance_left, ut_settings const& settings = ut_settings(),
        double distance_right = same_distance)
{
    test_one<Geometry, GeometryOut>(caseid, wkt, join_strategy, end_strategy,
        -1 ,-1, expected_area,
        distance_left, settings, distance_right);
}

template
<
    typename Geometry,
    typename GeometryOut,
    typename JoinStrategy,
    typename EndStrategy,
    typename DistanceStrategy,
    typename SideStrategy,
    typename PointStrategy
>
void test_with_custom_strategies(std::string const& caseid,
        std::string const& wkt,
        JoinStrategy const& join_strategy,
        EndStrategy const& end_strategy,
        DistanceStrategy const& distance_strategy,
        SideStrategy const& side_strategy,
        PointStrategy const& point_strategy,
        expectation_limits const& expected_area,
        ut_settings const& settings = ut_settings())
{
    namespace bg = boost::geometry;
    Geometry g;
    bg::read_wkt(wkt, g);
    bg::correct(g);

    bg::model::multi_polygon<GeometryOut> buffered;

    test_buffer<GeometryOut>
            (caseid, buffered, g,
            join_strategy, end_strategy,
            distance_strategy, side_strategy, point_strategy,
            expected_area, settings);
}

#endif
