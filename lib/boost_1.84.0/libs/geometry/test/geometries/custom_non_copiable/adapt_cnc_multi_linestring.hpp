// Boost.Geometry

// Copyright (c) 2023 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#ifndef GEOMETRY_TEST_ADAPT_CNC_MULTI_LINESTRING_HPP
#define GEOMETRY_TEST_ADAPT_CNC_MULTI_LINESTRING_HPP

#include "cnc_multi_linestring.hpp"

// 1. Adapt to Boost.Geometry
namespace boost { namespace geometry { namespace traits
{

template <typename Linestring>
struct tag<cnc_multi_linestring<Linestring>>
{
    using type = multi_linestring_tag;
};

template <typename Linestring>
struct clear<cnc_multi_linestring<Linestring> >
{
    static void apply(cnc_multi_linestring<Linestring>& multi)
    {
        multi.custom_clear();
    }
};

template <typename Linestring>
struct push_back<cnc_multi_linestring<Linestring> >
{
    static void apply(cnc_multi_linestring<Linestring>& multi, Linestring const& item)
    {
        multi.custom_push_back(item);
    }
    static inline void apply(cnc_multi_linestring<Linestring>& multi, Linestring&& item)
    {
        multi.custom_push_back_move(std::move(item));
    }
};

template <typename Linestring>
struct resize<cnc_multi_linestring<Linestring> >
{
    static void apply(cnc_multi_linestring<Linestring>& multi, std::size_t new_size)
    {
        multi.custom_resize(new_size);
    }
};


}}}

// 2a. Adapt to Boost.Range, meta-functions
namespace boost
{

template<typename Linestring>
struct range_mutable_iterator<cnc_multi_linestring<Linestring> >
{
    using type = decltype(std::declval<cnc_multi_linestring<Linestring>>().custom_begin());
};

template<typename Linestring>
struct range_const_iterator<cnc_multi_linestring<Linestring> >
{
    using type = decltype(std::declval<cnc_multi_linestring<Linestring>>().const_custom_begin());
};

}

// 2b. Adapt to Boost.Range, part 2, ADP
template<typename Linestring>
auto range_begin(cnc_multi_linestring<Linestring>& geo) { return geo.custom_begin(); }

template<typename Linestring>
auto range_begin(cnc_multi_linestring<Linestring> const& geo) { return geo.const_custom_begin(); }

template<typename Linestring>
auto range_end(cnc_multi_linestring<Linestring>& geo) { return geo.custom_end(); }

template<typename Linestring>
auto range_end(cnc_multi_linestring<Linestring> const& geo) { return geo.const_custom_end(); }

template<typename Linestring>
std::size_t range_calculate_size(cnc_multi_linestring<Linestring> const& geo) { return geo.custom_size(); }

#endif // GEOMETRY_TEST_ADAPT_CNC_MULTI_LINESTRING_HPP
