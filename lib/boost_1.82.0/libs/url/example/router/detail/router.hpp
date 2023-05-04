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
    BOOST_URL_DECL
    router_base();

    BOOST_URL_DECL
    virtual ~router_base();

    BOOST_URL_DECL
    void
    insert_impl(
        string_view s,
        any_resource const* v);

    BOOST_URL_DECL
    any_resource const*
    find_impl(
        segments_encoded_view path,
        string_view*& matches,
        string_view*& names) const noexcept;
};

} // detail
} // urls
} // boost

#endif
