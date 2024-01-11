// Boost.Geometry

// Copyright (c) 2023 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#ifndef GEOMETRY_TEST_ADAPT_CNC_POLYGON_HPP
#define GEOMETRY_TEST_ADAPT_CNC_POLYGON_HPP

#include "cnc_container.hpp"
#include "cnc_ring.hpp"
#include "cnc_polygon.hpp"

// 1. Adapt to Boost.Geometry
namespace boost { namespace geometry { namespace traits
{

template <typename Point>
struct tag<cnc_polygon<Point> >
{
    using type = polygon_tag;
};

template <typename Point>
struct ring_const_type<cnc_polygon<Point> >
{
    using type = typename cnc_polygon<Point>::custom_ring_type const&;
};

template <typename Point>
struct ring_mutable_type<cnc_polygon<Point> >
{
    using type = typename cnc_polygon<Point>::custom_ring_type&;
};

template <typename Point>
struct interior_const_type<cnc_polygon<Point> >
{
    using type = typename cnc_polygon<Point>::custom_int_type const& ;
};

template <typename Point>
struct interior_mutable_type<cnc_polygon<Point> >
{
    using type = typename cnc_polygon<Point>::custom_int_type&;
};


template <typename Point>
struct exterior_ring<cnc_polygon<Point> >
{
    static auto& get(cnc_polygon<Point>& p)
    {
        return p.custom_ext();
    }

    static auto const& get(cnc_polygon<Point> const& p)
    {
        return p.custom_ext();
    }
};

template <typename Point>
struct interior_rings<cnc_polygon<Point> >
{
    static auto& get(cnc_polygon<Point>& p)
    {
        return p.custom_int();
    }

    static auto const& get(cnc_polygon<Point> const& p)
    {
        return p.custom_int();
    }
};

}}} // namespace boost::geometry::traits

#endif // GEOMETRY_TEST_ADAPT_CNC_POLYGON_HPP
