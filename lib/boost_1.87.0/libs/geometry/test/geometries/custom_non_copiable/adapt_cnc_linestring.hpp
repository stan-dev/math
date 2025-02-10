// Boost.Geometry

// Copyright (c) 2023 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#ifndef GEOMETRY_TEST_ADAPT_CNC_LINESTRING_HPP
#define GEOMETRY_TEST_ADAPT_CNC_LINESTRING_HPP

#include "cnc_linestring.hpp"

// 1. Adapt to Boost.Geometry
namespace boost { namespace geometry { namespace traits
{

template <typename Point>
struct tag<cnc_linestring<Point>>
{
    using type = linestring_tag;
};

// Implement traits for mutable actions
// These are all optional traits (normally / default they are implemented
// conforming std:: functionality)

template <typename Point>
struct clear<cnc_linestring<Point> >
{
    static void apply(cnc_linestring<Point>& geo)
    {
        geo.custom_clear();
    }
};

template <typename Point>
struct push_back<cnc_linestring<Point> >
{
    static void apply(cnc_linestring<Point>& geo, Point const& point)
    {
        geo.custom_push_back(point);
    }
};

template <typename Point>
struct resize<cnc_linestring<Point> >
{
    static void apply(cnc_linestring<Point>& geo, std::size_t new_size)
    {
        geo.custom_resize(new_size);
    }
};


}}}

// 2a. Adapt to Boost.Range, meta-functions
namespace boost
{

template<typename Point>
struct range_mutable_iterator<cnc_linestring<Point> >
{
    using type = decltype(std::declval<cnc_linestring<Point>>().custom_begin());
};

template<typename Point>
struct range_const_iterator<cnc_linestring<Point> >
{
    using type = decltype(std::declval<cnc_linestring<Point>>().const_custom_begin());
};

}

// 2b. Adapt to Boost.Range, part 2, ADP
template<typename Point>
auto range_begin(cnc_linestring<Point>& geo) { return geo.custom_begin(); }

template<typename Point>
auto range_begin(cnc_linestring<Point> const& geo) { return geo.const_custom_begin(); }

template<typename Point>
auto range_end(cnc_linestring<Point>& geo) { return geo.custom_end(); }

template<typename Point>
auto range_end(cnc_linestring<Point> const& geo) { return geo.const_custom_end(); }

template<typename Point>
std::size_t range_calculate_size(cnc_linestring<Point> const& geo) { return geo.custom_size(); }

#endif // GEOMETRY_TEST_ADAPT_CNC_LINESTRING_HPP
