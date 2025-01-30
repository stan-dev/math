// Boost.Geometry

// Copyright (c) 2023 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#ifndef GEOMETRY_TEST_ADAPT_CNC_CONTAINER_HPP
#define GEOMETRY_TEST_ADAPT_CNC_CONTAINER_HPP

#include "cnc_container.hpp"

// 1. Adapt to Boost.Geometry
namespace boost { namespace geometry { namespace traits
{

template <typename Item>
struct clear<cnc_container<Item> >
{
    static void apply(cnc_container<Item>& container)
    {
        container.custom_clear();
    }
};

template <typename Item>
struct push_back<cnc_container<Item> >
{
    static void apply(cnc_container<Item>& container, Item const& item)
    {
        container.custom_push_back(item);
    }
    static void apply(cnc_container<Item>& container, Item&& item)
    {
        container.custom_push_back(std::move(item));
    }
};

template <typename Item>
struct resize<cnc_container<Item> >
{
    static void apply(cnc_container<Item>& container, std::size_t new_size)
    {
        container.custom_resize(new_size);
    }
};

}}} // namespace boost::geometry::traits


// 2a. Adapt to Boost.Range, meta-functions
namespace boost
{
    template<typename Item>
    struct range_mutable_iterator<cnc_container<Item> >
    {
        using type = decltype(std::declval<cnc_container<Item>>().custom_begin());
    };

    template<typename Item>
    struct range_const_iterator<cnc_container<Item> >
    {
        using type = decltype(std::declval<cnc_container<Item>>().const_custom_begin());
    };


} // namespace boost


// 2b. Adapt to Boost.Range, part 2, ADP
template<typename Item>
auto range_begin(cnc_container<Item>& container) { return container.custom_begin(); }

template<typename Item>
auto range_begin(cnc_container<Item> const& container) { return container.const_custom_begin(); }

template<typename Item>
auto range_end(cnc_container<Item>& container) { return container.custom_end(); }

template<typename Item>
auto range_end(cnc_container<Item> const& container) { return container.const_custom_end(); }

template<typename Item>
std::size_t range_calculate_size(cnc_container<Item> const& container) { return container.custom_size(); }





#endif // GEOMETRY_TEST_ADAPT_CNC_CONTAINER_HPP
