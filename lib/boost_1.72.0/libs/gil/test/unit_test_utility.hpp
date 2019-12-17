//
// Copyright 2019 Mateusz Loskot <mateusz at loskot dot net>
//
// Distribtted under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//
#ifndef BOOST_GIL_TEST_UNIT_TEST_UTILITY_HPP
#define BOOST_GIL_TEST_UNIT_TEST_UTILITY_HPP

#include <boost/gil/color_base_algorithm.hpp> // static_for_each
#include <boost/gil/pixel.hpp>
#include <boost/gil/promote_integral.hpp>

#include <boost/core/typeinfo.hpp>

#include <cstdint>
#include <ostream>

namespace boost { namespace gil {

namespace detail {

struct print_color_base
{
    std::ostream& os_;
    std::size_t element_index_{0};
    print_color_base(std::ostream& os) : os_(os) {}

    template <typename Element>
    void operator()(Element& c)
    {
        typename ::boost::gil::promote_integral<Element>::type const v(c);
        if (element_index_ > 0) os_ << ", ";
        os_ << "v" << element_index_ << "=" << v;
        ++element_index_;
    }
};

} // namespace detail

// Make `point` printable for BOOST_TEST()
template <typename T>
std::ostream& operator<<(std::ostream& os, point<T> const& p)
{
    os << "point<" << boost::core::demangled_name(typeid(T)) << ">";
    os << "(" << p.x << ", " << p.y << ")" << std::endl;
    return os;
}

// Make `pixel` printable for BOOST_TEST()
template <typename ChannelValue, typename Layout>
std::ostream& operator<<(std::ostream& os, pixel<ChannelValue, Layout> const& p)
{
    os << "pixel<"
       << "Channel=" << boost::core::demangled_name(typeid(ChannelValue))
       << ", Layout=" << boost::core::demangled_name(typeid(Layout))
       << ">(";

    static_for_each(p, detail::print_color_base{os});
    os << ")" << std::endl;
    return os;
}

}} // namespace boost::gil

#endif
