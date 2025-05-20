//
// Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

#ifndef BOOST_URL_DETAIL_ROUTER_HPP
#define BOOST_URL_DETAIL_ROUTER_HPP

#include <boost/url/pct_string_view.hpp>
#include <boost/url/segments_encoded_view.hpp>
#include <boost/url/grammar/delim_rule.hpp>
#include <boost/url/grammar/optional_rule.hpp>
#include <boost/url/grammar/range_rule.hpp>
#include <boost/url/grammar/tuple_rule.hpp>
#include <string>

namespace boost {
namespace urls {
namespace detail {

class router_base
{
    void* impl_{nullptr};

public:
    // A type-erased router resource
    struct any_resource
    {
        virtual ~any_resource() = default;
        virtual void const* get() const noexcept = 0;
    };

protected:
    router_base();

    virtual ~router_base();

    void
    insert_impl(
        core::string_view s,
        any_resource const* v);

    any_resource const*
    find_impl(
        segments_encoded_view path,
        core::string_view*& matches,
        core::string_view*& names) const noexcept;
};

} // detail
} // urls
} // boost

#endif
