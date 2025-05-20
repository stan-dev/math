// Boost.Geometry

// Copyright (c) 2023 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#ifndef GEOMETRY_TEST_CNC_POLYGON_HPP
#define GEOMETRY_TEST_CNC_POLYGON_HPP

#include "cnc_container.hpp"
#include "cnc_ring.hpp"

template <typename P>
class cnc_polygon
{
public :
    using custom_ring_type = cnc_ring<P> ;
    using custom_int_type = cnc_container<custom_ring_type>;

    auto& custom_ext() { return m_ext; }
    auto& custom_int() { return m_int; }

    auto const& custom_ext() const { return m_ext; }
    auto const& custom_int() const { return m_int; }

    // Make sure the polygon can only be moved and never copied.
    cnc_polygon() = default;
    cnc_polygon(cnc_polygon&& g) = default;
    ~cnc_polygon() = default;
    cnc_polygon(const cnc_polygon& g) = delete;
    cnc_polygon& operator=(const cnc_polygon&) = delete;

private :
    custom_ring_type m_ext;
    custom_int_type m_int;
};

#endif // GEOMETRY_TEST_CNC_POLYGON_HPP
