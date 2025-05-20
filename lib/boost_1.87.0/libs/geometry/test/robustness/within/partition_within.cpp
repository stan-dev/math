// Boost.Geometry
//
// Copyright (c) 2023 Barend Gehrels, Amsterdam, the Netherlands.
//
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/geometry.hpp>

#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/register/point.hpp>

#if defined(TEST_WITH_SVG)
# include <boost/geometry/io/svg/svg_mapper.hpp>
#endif

#include <chrono>
#include <random>
#include <fstream>

namespace bg = boost::geometry;

struct point_item
{
    std::size_t id = 0;
    double x;
    double y;
};

template <typename Ring>
struct ring_item
{
    using point_t = typename bg::point_type<Ring>::type;
    using ring_t = Ring;
    std::size_t id = 0;
    bg::model::box<point_t> box;
    bg::model::ring<point_t> ring;
};


BOOST_GEOMETRY_REGISTER_POINT_2D(point_item, double, cs::cartesian, x, y)


struct expand_for_point
{
    template <typename Box, typename InputItem>
    static inline void apply(Box& total, InputItem const& item)
    {
        bg::expand(total, item);
    }
};

struct overlaps_point
{
    template <typename Box, typename InputItem>
    static inline bool apply(Box const& box, InputItem const& item)
    {
        return ! bg::disjoint(item, box);
    }
};


struct expand_for_ring
{
    template <typename Box, typename InputItem>
    static inline void apply(Box& total, InputItem const& item)
    {
        bg::expand(total, item.box);
    }
};

struct overlaps_ring
{
    template <typename Box, typename InputItem>
    static inline bool apply(Box const& box, InputItem const& item)
    {
        typename bg::strategy::disjoint::services::default_strategy
            <
                Box, Box
            >::type strategy;

        return ! bg::detail::disjoint::disjoint_box_box(box, item.box, strategy);
    }
};

struct point_in_ring_visitor
{
    std::size_t count = 0;

    template <typename Point, typename BoxItem>
    inline bool apply(Point const& point, BoxItem const& ring_item)
    {
        if (bg::within(point, ring_item.ring))
        {
            count++;
        }
        return true;
    }
};

#if defined(TEST_WITH_SVG)
template <typename Box, typename Points, typename Rings>
struct svg_visitor
{
    std::vector<Box> boxes;
    Points const& m_points;
    Rings const& m_rings;
    std::size_t m_size = 0;
    std::size_t m_index = 0;

    svg_visitor(std::size_t size, Points const& points, Rings const& rings)
        : m_points(points)
        , m_rings(rings)
        , m_size(size)
    {}

    inline void apply(Box const& box, int level)
    {
        std::ostringstream filename;
        filename << "partition_demo_" << std::setfill('0') << std::setw(3) << m_index++ << "_" << level << ".svg";
        std::ofstream svg(filename.str());

        bg::svg_mapper<point_item> mapper(svg, 800, 800);

        {
            point_item p;
            p.x = -1; p.y = -1; mapper.add(p);
            p.x = m_size + 1; p.y = m_size + 1; mapper.add(p);
        }

        for (auto const& item : m_rings)
        {
            mapper.map(item.ring, "opacity:0.6;fill:rgb(0,255,0);stroke:rgb(0,0,0);stroke-width:0.1");
        }
        for (auto const& point : m_points)
        {
            mapper.map(point, "fill:rgb(0,0,255);stroke:rgb(0,0,100);stroke-width:0.1", 3);
        }

        for (auto const& b : boxes)
        {
            mapper.map(b, "fill:none;stroke-width:2;stroke:rgb(64,64,64);");
        }
        mapper.map(box, "fill:none;stroke-width:4;stroke:rgb(255, 0, 0);");

        boxes.push_back(box);
    }
};
#endif


template <typename Collection>
void fill_points(Collection& collection, std::size_t size, std::size_t count)
{
    std::random_device rd;

    std::default_random_engine rde(rd());
    std::uniform_real_distribution<double> uniform_dist(0, size - 1);
    int const mean_x = uniform_dist(rde);
    int const mean_y = uniform_dist(rde);
 
    // Generate a normal distribution around these means
    std::seed_seq seed2{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    std::mt19937 e2(seed2);
    std::normal_distribution<> normal_dist_x(mean_x, size / 5.0);
    std::normal_distribution<> normal_dist_y(mean_y, size / 5.0);

    int n = 0;
    for (int i = 0; n < count && i < count * count; i++)
    {
        double const x = normal_dist_x(e2);
        double const y = normal_dist_y(e2);
        if (x >= 0 && y >= 0 && x < size && y < size)
        {
            typename boost::range_value<Collection>::type item;
            item.x = x;
            item.y = y;
            collection.push_back(item);
            n++;
        }
    }
}

template <typename Ring>
auto create_ring(std::default_random_engine& engine)
{
    // For the initial (unmoved, unsized) polygons, in [0..1]
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    auto cross_product = [](auto p1, auto p2, auto p3) { return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x); };

    auto is_concave = [&cross_product](auto const& points) {
        int const point_count = points.size();
        if (point_count < 3)
        {
            return false;
        }
        for (int i = 0; i < point_count; i++)
        {
            if (cross_product(points[i], points[(i + 1) % point_count], points[(i + 2) % point_count]) > 0.0)
            {
                return true;
            }
        }
        return false;
    };

    std::size_t iteration = 0;
    Ring ring;
    while (ring.size() < 5 && iteration < 100)
    {
        double const xp = distribution(engine);
        double const yp = distribution(engine);
        point_item const p{0, xp, yp};
        ring.push_back(p);
        if (is_concave(ring))
        {
            ring.pop_back();
        }
        iteration++;
    }
    // Close it and make it clockwise
    bg::correct(ring);
    return ring;
}


template <typename Collection>
void fill_rings(Collection& collection, std::size_t size, std::size_t count)
{
    using item_t = typename boost::range_value<Collection>::type;
    using ring_t = typename item_t::ring_t;

    // For the size of the polygons (w/h)
    double const min_dimension = 5.0;
    double const max_dimension = size / 15.0;
    if (max_dimension <= min_dimension)
    {
        throw std::runtime_error("Size is too small");
    }

    std::random_device rd;
    std::default_random_engine dre(rd());

    // For the polygon dimensions
    std::uniform_real_distribution<double> uniform_dist_dimension(min_dimension, max_dimension);

    // For the polygon location
    std::uniform_real_distribution<double> uniform_dist_location(0.0, size - min_dimension);


    int n = 0;
    for (int i = 0; n < count && i < count * count; i++)
    {
        // Generate polygon location (x,y) and dimension (w,h)
        double const w = uniform_dist_dimension(dre);
        double const h = uniform_dist_dimension(dre);

        double const x = uniform_dist_location(dre);
        double const y = uniform_dist_location(dre);
        if (x + w >= size || y + h >= size)
        {
            continue;
        }

        item_t item;
        item.id = n + 1;

        item.ring = create_ring<ring_t>(dre);

        // Avoid small oblong slivers by having a minimum size.
        if (bg::area(item.ring) > 0.2)
        {
            // Increase the polygon size
            bg::for_each_point(item.ring, [&w, &h](auto& point) { point.x *= w; point.y *= h;} );

            // Move the polygon to (x,y)
            bg::for_each_point(item.ring, [&x, &y](auto& point) { point.x += x; point.y += y;} );

            // Calculate its box
            bg::envelope(item.ring, item.box);

            collection.push_back(item);
            n++;
        }
    }
}

void call_within(std::size_t size, std::size_t count)
{
    using box_type = bg::model::box<point_item>;
    using ring_type = bg::model::ring<point_item>;
    std::vector<point_item> points;
    std::vector<ring_item<ring_type>> rings;

    fill_points(points, size, count);
    fill_rings(rings, size, count);

    auto report = [&points, &rings](const char* title, auto const& start, std::size_t within_count)
    {
        auto const finish = std::chrono::steady_clock::now();
        double const elapsed 
            = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start).count();
        std::cout << title << " time: " << std::setprecision(6)
                << elapsed / 1000000.0 << " ms" << std::endl;
        std::cout << "Points in rings: " << within_count
                << " of " << points.size() << " / " << rings.size() << std::endl;
        return elapsed;
    };

    point_in_ring_visitor count_visitor;
    {
#if defined(TEST_WITH_SVG)

        using partition_box_visitor_type = svg_visitor<box_type, std::vector<point_item>, std::vector<ring_item<ring_type>>>;
        partition_box_visitor_type partition_box_visitor(size, points, rings);
#else
        using partition_box_visitor_type = bg::detail::partition::visit_no_policy;
        partition_box_visitor_type partition_box_visitor;
#endif

        auto const start = std::chrono::steady_clock::now();

        bg::partition
            <
                box_type,
                bg::detail::partition::include_all_policy,
                bg::detail::partition::include_all_policy
            >::apply(points, rings, count_visitor, expand_for_point(), overlaps_point(),
                    expand_for_ring(), overlaps_ring(), 16, partition_box_visitor);
        report("Partition", start, count_visitor.count);
    }

    {
        // Verify and compare it with a quadratic loop
        auto const start = std::chrono::steady_clock::now();
        std::size_t count = 0;
        for (auto const& point : points)
        {
            for (auto const& ring : rings)
            {
                if (bg::within(point, ring.ring))
                {
                    count++;
                }
            }
        }
        report("Quadratic loop", start, count);

        if (count != count_visitor.count)
        {
            std::cerr << "ERROR: counts are not equal" << std::endl;
        }
    }
}

int main()
{
    for (int i = 0; i < 10; i++)
    {
        call_within(100, 2000);
    }
    call_within(200, 20000);

    return 0;
}
