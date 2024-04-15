//
// Copyright (c) 2023 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

#ifndef BOOST_URL_RFC_DETAIL_IMPL_IPV6_ADDRZ_RULE_IPP
#define BOOST_URL_RFC_DETAIL_IMPL_IPV6_ADDRZ_RULE_IPP

#include <boost/url/detail/config.hpp>
#include <boost/url/grammar/parse.hpp>
#include <boost/url/rfc/detail/ipv6_addrz_rule.hpp>
#include <boost/url/rfc/ipv6_address_rule.hpp>
#include <boost/url/rfc/unreserved_chars.hpp>
#include <boost/url/rfc/pct_encoded_rule.hpp>

namespace boost {
namespace urls {
namespace detail {

auto
ipv6_addrz_rule_t::
parse(
    char const*& it,
    char const* const end
        ) const noexcept ->
    system::result<value_type>
{
    value_type t;
    auto rv1 = grammar::parse(
        it, end, ipv6_address_rule);
    if (! rv1)
        return rv1.error();
    t.ipv6 = *rv1;

    // "%25"
    auto it0 = it;
    if (end - it < 3 ||
        *it != '%' ||
        *(it + 1) != '2' ||
        *(it + 2) != '5')
    {
        BOOST_URL_RETURN_EC(
            grammar::error::invalid);
    }
    it += 3;

    // ZoneID = 1*( unreserved / pct-encoded )
    // Parse as many (unreserved / pct-encoded)
    // as available
    auto rv2 = grammar::parse(
            it, end,
            pct_encoded_rule(unreserved_chars));
    if(!rv2 || rv2->empty())
    {
        it = it0;
        BOOST_URL_RETURN_EC(
            grammar::error::invalid);
    }
    else
    {
        t.zone_id = *rv2;
    }
    return t;
}

} // detail
} // urls
} // boost

#endif
