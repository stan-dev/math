//
// Copyright (c) 2016-2019 Vinnie Falco (vinnie dot falco at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//


#include <boost/url/detail/config.hpp>
#include <boost/url/rfc/uri_rule.hpp>
#include <boost/url/rfc/query_rule.hpp>
#include "detail/fragment_part_rule.hpp"
#include "detail/hier_part_rule.hpp"
#include "detail/query_part_rule.hpp"
#include "detail/scheme_rule.hpp"
#include <boost/url/grammar/delim_rule.hpp>
#include <boost/url/grammar/tuple_rule.hpp>
#include <boost/url/grammar/optional_rule.hpp>
#include <boost/url/grammar/parse.hpp>

namespace boost {
namespace urls {

auto
implementation_defined::uri_rule_t::
parse(
    char const*& it,
    char const* const end
        ) const noexcept ->
    system::result<value_type>
{
    detail::url_impl u(detail::url_impl::from::string);
    u.cs_ = it;

    // scheme
    {
        auto rv = grammar::parse(
            it, end,
            grammar::tuple_rule(
                detail::scheme_rule(),
                grammar::squelch(
                    grammar::delim_rule(':'))));
        if(! rv)
            return rv.error();
        u.apply_scheme(rv->scheme);
    }

    // hier_part
    {
        auto rv = grammar::parse(
            it, end,
            detail::hier_part_rule);
        if(! rv)
            return rv.error();
        if(rv->has_authority)
            u.apply_authority(rv->authority);
        u.apply_path(
            rv->path,
            rv->segment_count);
    }

    // [ "?" query ]
    {
        auto rv = grammar::parse(
            it, end, detail::query_part_rule);
        if(! rv)
            return rv.error();
        if(rv->has_query)
        {
            // map "?" to { {} }
            u.apply_query(
                rv->query,
                rv->count);
        }
    }

    // [ "#" fragment ]
    {
        auto rv = grammar::parse(
            it, end, detail::fragment_part_rule);
        if(! rv)
            return rv.error();
        if(rv->has_fragment)
            u.apply_frag(rv->fragment);
    }

    return u.construct();
}

} // urls
} // boost

