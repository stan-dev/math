// Boost.Geometry

// Copyright (c) 2023 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#ifndef GEOMETRY_TEST_CNC_CONTAINER_HPP
#define GEOMETRY_TEST_CNC_CONTAINER_HPP

#include <cstddef>
#include <deque>

// Define a custom class which is only movable and not copiable, and cannot be indexed.
// Its derivatives will be linestring, ring, multi*
// For readability purposes we use the "cnc" abbreviation as "custom_non_copiable"
template <typename Item>
class cnc_container
{
public :

    // Make sure the library does not use, for example, begin/end or indexing
    // Therefor all functions are custom functions.
    std::size_t custom_size() const { return m_private_deque.size(); }

    auto const_custom_begin() const { return m_private_deque.begin(); }
    auto const_custom_end() const { return m_private_deque.end(); }
    auto custom_begin() { return m_private_deque.begin(); }
    auto custom_end() { return m_private_deque.end(); }

    void custom_clear() { m_private_deque.clear(); }
    void custom_push_back(Item const& p) { m_private_deque.push_back(p); }
    void custom_push_back_move(Item&& p) { m_private_deque.push_back(std::move(p)); }
    void custom_resize(std::size_t new_size) { m_private_deque.resize(new_size); }

    auto& custom_get(std::size_t index) { return m_private_deque[index]; }

    // Make sure the object, and inherited objects, can only be moved and never copied.
    cnc_container() = default;
    cnc_container(cnc_container&& g) = default;
    ~cnc_container() = default;
    cnc_container(const cnc_container& g) = delete;
    cnc_container& operator=(const cnc_container&) = delete;

private :
    std::deque<Item> m_private_deque;
};

#endif // GEOMETRY_TEST_CNC_CONTAINER_HPP
