// Boost.Geometry

// Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.

// This file was modified by Oracle on 2019-2024.
// Modifications copyright (c) 2019-2024 Oracle and/or its affiliates.
// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <geometry_test_common.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/optional.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/register/point.hpp>

namespace bg = boost::geometry;

#if defined(TEST_WITH_SVG)
#include "test_buffer_svg.hpp"
#endif

// A point with extra info, such that
// - it can influence the buffer with dynamically (input)
// - it can receive the side at which the buffer was (output)
struct specific_point
{
    double x{0.0};
    double y{0.0};
    double distance_left{0.0};
    double distance_right{0.0};
    boost::optional<bg::strategy::buffer::buffer_side_selector> side{};
};

BOOST_GEOMETRY_REGISTER_POINT_2D(specific_point, double, cs::cartesian, x, y)

struct specific_buffer_side_strategy
{
    bool equidistant() const
    {
        return true;
    }

    template <typename Point, typename OutputRange, typename DistanceStrategy>
    bg::strategy::buffer::result_code apply(Point const& input_p1, Point const& input_p2,
                                            bg::strategy::buffer::buffer_side_selector side,
                                            DistanceStrategy const& distance,
                                            OutputRange &output_range) const
    {
        const auto result =
            bg::strategy::buffer::side_straight::apply(input_p1, input_p2, side, distance, output_range);

        for (auto& point : output_range)
        {
            point.side = side;
        }
        return result;
    }
};

struct specific_buffer_distance_strategy
{
    /// Returns the distance-value.
    template <typename Point>
    double apply(Point const& p, Point const& q,
                 bg::strategy::buffer::buffer_side_selector side) const
    {
        return bg::strategy::buffer::buffer_side_left == side
                   ? std::max(p.distance_left, q.distance_left)
                   : std::max(p.distance_right, q.distance_right);
    }

    bool negative() const
    {
        return false;
    }

    bool empty(bg::strategy::buffer::buffer_side_selector) const
    {
        return false;
    }

    /// Returns the max distance distance up to the buffer will reach.
    template <typename JoinStrategy, typename EndStrategy>
    double max_distance(JoinStrategy const &, EndStrategy const &) const
    {
        constexpr double unit_buffer_width{3.0};
        return 6.0 * unit_buffer_width;
    }

    /// Returns the distance at which the input is simplified before the buffer process.
    double simplify_distance() const
    {
        constexpr double unit_buffer_width{3.0};
        return unit_buffer_width / 1000.0;
    }
};

void test_buffer(std::string const& caseid, std::string const& wkt, std::vector<double> const& widths, double expected_area)
{
    using ls_t = boost::geometry::model::linestring<specific_point>;
    using mls_t = boost::geometry::model::multi_linestring<ls_t>;
    using polygon_t = boost::geometry::model::polygon<specific_point>;

    mls_t mls;
    boost::geometry::read_wkt(wkt, mls);

    if (mls.size() != widths.size())
    {
        throw std::runtime_error("There should be correct widths");
    }

    using point_type = specific_point;
    using strategy_t = bg::strategies::buffer::services::default_strategy<
        mls_t>::type;
    strategy_t strategy;

#if defined(TEST_WITH_SVG)
    bg::model::box<point_type> envelope;
    bg::envelope(mls, envelope, strategy);

    buffer_svg_mapper<point_type> buffer_mapper(caseid);

    std::ostringstream filename;
    filename << "/tmp/buffer_variable_width_" << caseid << ".svg";
    std::ofstream svg(filename.str().c_str());
    typedef bg::svg_mapper<point_type> mapper_type;
    mapper_type mapper(svg, 1000, 800);

    svg_visitor<mapper_type, bg::model::box<point_type>> visitor(mapper);

    // Set the SVG boundingbox, with a margin. The margin is necessary because
    // drawing is already started before the buffer is finished. It is not
    // possible to "add" the buffer (unless we buffer twice).
    buffer_mapper.prepare(mapper, visitor, envelope, 2.0);
#else
    bg::detail::buffer::visit_pieces_default_policy visitor;
#endif

    for (std::size_t i = 0; i < mls.size(); i++)
    {
        for (auto& point : mls[i])
        {
            point.distance_left = widths[i];
            point.distance_right = widths[i];
        }
    }

    specific_buffer_side_strategy buffer_side_strategy;
    specific_buffer_distance_strategy distance_strategy;
    boost::geometry::strategy::buffer::join_miter join_strategy;
    boost::geometry::strategy::buffer::end_flat end_strategy;
    boost::geometry::strategy::buffer::point_circle point_strategy;
    boost::geometry::model::multi_polygon<polygon_t> result;
    boost::geometry::detail::buffer::buffer_inserter
        <
            polygon_t
        >(mls,
          std::back_inserter(result),
          distance_strategy, buffer_side_strategy, join_strategy,
          end_strategy, point_strategy,
          strategy,
          visitor);

#if defined(TEST_WITH_SVG)
    buffer_mapper.map_input_output(mapper, mls, result, false);
#endif

    std::cout << caseid << " : " << boost::geometry::area(result) << std::endl;
    BOOST_CHECK_CLOSE_FRACTION(boost::geometry::area(result), expected_area, 0.001);
}

int test_main(int argc, char **argv)
{
    test_buffer("case1", "MULTILINESTRING((-2.99999999998012700785920969792641699314117431640625 10.4999999997948929575386500800959765911102294921875,18.00000000002045652536253328435122966766357421875 10.4999999997948947338954894803464412689208984375),(-2.99999999998012700785920969792641699314117431640625 10.4999999997948929575386500800959765911102294921875,17.99999999995061017443731543608009815216064453125 -4.5000000001943778471513724070973694324493408203125),(-2.99999999998012700785920969792641699314117431640625 10.4999999997948929575386500800959765911102294921875,11.99999999996175148453403380699455738067626953125 -10.49999999983214848953139153309166431427001953125))",
             {4.5, 1.5, 1.5}, 304.294);
    test_buffer("case2", "MULTILINESTRING((-18.00000000002045652536253328435122966766357421875 10.4999999997948929575386500800959765911102294921875,-2.99999999998012700785920969792641699314117431640625 10.4999999997948929575386500800959765911102294921875),(-2.99999999998012700785920969792641699314117431640625 10.4999999997948929575386500800959765911102294921875,17.99999999995061017443731543608009815216064453125 -4.5000000001943778471513724070973694324493408203125),(-2.99999999998012700785920969792641699314117431640625 10.4999999997948929575386500800959765911102294921875,11.99999999996175148453403380699455738067626953125 -10.49999999983214848953139153309166431427001953125))",
             {6.0, 1.5, 1.5}, 319.767);
    test_buffer("case3", "MULTILINESTRING((-18.00000000002045652536253328435122966766357421875 10.4999999997948929575386500800959765911102294921875,-2.99999999998012700785920969792641699314117431640625 10.4999999997948929575386500800959765911102294921875),(-2.99999999998012700785920969792641699314117431640625 10.4999999997948929575386500800959765911102294921875,18.00000000002045652536253328435122966766357421875 10.4999999997948947338954894803464412689208984375),(-2.99999999998012700785920969792641699314117431640625 10.4999999997948929575386500800959765911102294921875,11.99999999996175148453403380699455738067626953125 -10.49999999983214848953139153309166431427001953125))",
             {6.0, 4.5, 1.5}, 429.831);
    test_buffer("case4", "MULTILINESTRING((-18.00000000002045652536253328435122966766357421875 10.4999999997948929575386500800959765911102294921875,-2.99999999998012700785920969792641699314117431640625 10.4999999997948929575386500800959765911102294921875),(-2.99999999998012700785920969792641699314117431640625 10.4999999997948929575386500800959765911102294921875,18.00000000002045652536253328435122966766357421875 10.4999999997948947338954894803464412689208984375),(-2.99999999998012700785920969792641699314117431640625 10.4999999997948929575386500800959765911102294921875,17.99999999995061017443731543608009815216064453125 -4.5000000001943778471513724070973694324493408203125))",
             {6.0, 4.5, 1.5}, 423.195);
    return 0;
}
