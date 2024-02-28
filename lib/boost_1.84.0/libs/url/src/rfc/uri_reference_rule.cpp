//
// Copyright (c) 2016-2019 Vinnie Falco (vinnie dot falco at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

#ifndef BOOST_URL_RFC_IMPL_URI_REFERENCE_RULE_IPP
#define BOOST_URL_RFC_IMPL_URI_REFERENCE_RULE_IPP

#include <boost/url/detail/config.hpp>
#include <boost/url/rfc/uri_reference_rule.hpp>
#include <boost/url/rfc/uri_rule.hpp>
#include <boost/url/rfc/relative_ref_rule.hpp>
#include <boost/url/grammar/parse.hpp>
#include <boost/url/grammar/variant_rule.hpp>
#include <boost/variant2/variant.hpp>

namespace boost {
namespace urls {

auto
uri_reference_rule_t::
parse(
    char const*& it,
    char const* const end
        ) const noexcept ->
    system::result<value_type>
{
    auto rv = grammar::parse(
        it, end,
        grammar::variant_rule(
            uri_rule,
            relative_ref_rule));
    if(! rv)
        return rv.error();
    switch(rv->index())
    {
    default:
    case 0:
        return boost::variant2::get<0>(*rv);
    case 1:
        return boost::variant2::get<1>(*rv);
    }
}

} // urls
} // boost

#endif
