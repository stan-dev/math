// Boost.Geometry
// Unit Test Helper

// Copyright (c) 2022 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_GEOMETRY_TEST_ALGORITHMS_BUFFER_DEBUG_BUFFER_INFO_HPP
#define BOOST_GEOMETRY_TEST_ALGORITHMS_BUFFER_DEBUG_BUFFER_INFO_HPP

namespace boost { namespace geometry { namespace debug
{

// Helper method to translate enumeration to character
inline char piece_type_char(strategy::buffer::piece_type const& type)
{
    using namespace strategy::buffer;
    switch(type)
    {
        case buffered_segment : return 's';
        case buffered_join : return 'j';
        case buffered_round_end : return 'r';
        case buffered_flat_end : return 'f';
        case buffered_point : return 'p';
        case buffered_concave : return 'c';
        default : return '?';
    }
}

}}}

#endif
