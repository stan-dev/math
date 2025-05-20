// Boost.Geometry

// Copyright (c) 2023 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#ifndef GEOMETRY_TEST_ADAPT_CNC_MULTI_POLYGON_HPP
#define GEOMETRY_TEST_ADAPT_CNC_MULTI_POLYGON_HPP

#include "cnc_multi_polygon.hpp"

// 1. Adapt to Boost.Geometry
namespace boost { namespace geometry { namespace traits
{

template <typename Polygon>
struct tag<cnc_multi_polygon<Polygon>>
{
    using type = multi_polygon_tag;
};

template <typename Polygon>
struct clear<cnc_multi_polygon<Polygon> >
{
    static void apply(cnc_multi_polygon<Polygon>& multi)
    {
        multi.custom_clear();
    }
};

template <typename Polygon>
struct push_back<cnc_multi_polygon<Polygon> >
{
    static void apply(cnc_multi_polygon<Polygon>& multi, Polygon const& item)
    {
        multi.custom_push_back(item);
    }
    static inline void apply(cnc_multi_polygon<Polygon>& multi, Polygon&& item)
    {
        multi.custom_push_back_move(std::move(item));
    }
};

template <typename Polygon>
struct resize<cnc_multi_polygon<Polygon> >
{
    static void apply(cnc_multi_polygon<Polygon>& multi, std::size_t new_size)
    {
        multi.custom_resize(new_size);
    }
};


}}}

// 2a. Adapt to Boost.Range, meta-functions
namespace boost
{

template<typename Polygon>
struct range_mutable_iterator<cnc_multi_polygon<Polygon> >
{
    using type = decltype(std::declval<cnc_multi_polygon<Polygon>>().custom_begin());
};

template<typename Polygon>
struct range_const_iterator<cnc_multi_polygon<Polygon> >
{
    using type = decltype(std::declval<cnc_multi_polygon<Polygon>>().const_custom_begin());
};

}

// 2b. Adapt to Boost.Range, part 2, ADP
template<typename Polygon>
auto range_begin(cnc_multi_polygon<Polygon>& geo) { return geo.custom_begin(); }

template<typename Polygon>
auto range_begin(cnc_multi_polygon<Polygon> const& geo) { return geo.const_custom_begin(); }

template<typename Polygon>
auto range_end(cnc_multi_polygon<Polygon>& geo) { return geo.custom_end(); }

template<typename Polygon>
auto range_end(cnc_multi_polygon<Polygon> const& geo) { return geo.const_custom_end(); }

template<typename Polygon>
std::size_t range_calculate_size(cnc_multi_polygon<Polygon> const& geo) { return geo.custom_size(); }

#endif // GEOMETRY_TEST_ADAPT_CNC_MULTI_POLYGON_HPP
