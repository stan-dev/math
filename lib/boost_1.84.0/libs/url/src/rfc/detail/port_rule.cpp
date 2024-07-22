//
// Copyright (c) 2016-2019 Vinnie Falco (vinnie dot falco at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

#ifndef BOOST_URL_IMPL_PORT_RULE_IPP
#define BOOST_URL_IMPL_PORT_RULE_IPP

#include <boost/url/detail/config.hpp>
#include <boost/url/rfc/detail/port_rule.hpp>
#include <boost/url/grammar/parse.hpp>
#include <boost/url/grammar/token_rule.hpp>
#include <boost/url/grammar/unsigned_rule.hpp>
#include <boost/static_assert.hpp>
#include <type_traits>

namespace boost {
namespace urls {
namespace detail {

auto
port_rule::
parse(
    char const*& it,
    char const* end) const noexcept ->
        system::result<value_type>
{
    value_type t;
    auto const start = it;
    while(
        it != end &&
        *it == '0')
    {
        ++it;
    }

    if (it != end)
    {
        grammar::unsigned_rule<std::uint16_t> r;
        auto it0 = it;
        auto rv = r.parse(it, end);
        if (rv)
        {
            // number < max uint16_t
            t.str = core::string_view(start, it);
            t.has_number = true;
            t.number = *rv;
            return t;
        }
        it = it0;
        if (grammar::digit_chars(*it))
        {
            // number > max uint16_t
            while (
                it != end &&
                grammar::digit_chars(*it))
            {
                ++it;
            }
            t.str = core::string_view(start, it);
            t.has_number = true;
            t.number = 0;
            return t;
        }
    }
    // no digits
    t.str = core::string_view(start, it);
    t.has_number = it != end;
    t.number = 0;
    return t;
}

auto
port_part_rule_t::
parse(
    char const*& it,
    char const* end) const noexcept ->
        system::result<value_type>
{
    value_type t;
    if( it == end ||
        *it != ':')
    {
        t.has_port = false;
        return t;
    }
    ++it;
    auto rv = grammar::parse(
        it, end, port_rule{});
    if(! rv)
        return rv.error();
    t.has_port = true;
    t.port = rv->str;
    t.has_number = rv->has_number;
    t.port_number = rv->number;
    return t;
}

} // detail
} // urls
} // boost

#endif
