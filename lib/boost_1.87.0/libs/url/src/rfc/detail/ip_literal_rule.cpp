//
// Copyright (c) 2016-2019 Vinnie Falco (vinnie dot falco at gmail dot com)
// Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//


#include <boost/url/detail/config.hpp>
#include <boost/url/rfc/ipv6_address_rule.hpp>
#include "ip_literal_rule.hpp"
#include "ipv6_addrz_rule.hpp"
#include <boost/url/grammar/delim_rule.hpp>
#include <boost/url/grammar/parse.hpp>
#include <boost/url/grammar/tuple_rule.hpp>
#include "ipvfuture_rule.hpp"

namespace boost {
namespace urls {
namespace detail {

auto
ip_literal_rule_t::
parse(
    char const*& it,
    char const* const end
        ) const noexcept ->
    system::result<value_type>
{
    value_type t;

    // '['
    {
        auto rv = grammar::parse(
            it, end, grammar::delim_rule('['));
        if(! rv)
            return rv.error();
    }
    if(it == end)
    {
        // end
        BOOST_URL_RETURN_EC(
            grammar::error::invalid);
    }
    if(*it != 'v')
    {
        // IPv6address
        auto it0 = it;
        auto rv = grammar::parse(
            it, end,
            grammar::tuple_rule(
                ipv6_address_rule,
                grammar::squelch(
                    grammar::delim_rule(']'))));
        if(! rv)
        {
            // IPv6addrz
            it = it0;
            auto rv2 = grammar::parse(
                it, end,
                grammar::tuple_rule(
                    ipv6_addrz_rule,
                    grammar::squelch(
                        grammar::delim_rule(']'))));
            if(! rv2)
                return rv2.error();
            t.ipv6 = rv2->ipv6;
            t.is_ipv6 = true;
            return t;
        }
        t.ipv6 = *rv;
        t.is_ipv6 = true;
        return t;
    }
    {
        // IPvFuture
        auto rv = grammar::parse(
            it, end, 
            grammar::tuple_rule(
                ipvfuture_rule,
                grammar::squelch(
                    grammar::delim_rule(']'))));
        if(! rv)
            return rv.error();
        t.is_ipv6 = false;
        t.ipvfuture = rv->str;
        return t;
    }
}

} // detail
} // urls
} // boost

