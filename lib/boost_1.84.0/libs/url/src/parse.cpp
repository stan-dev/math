//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
// Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

#ifndef BOOST_URL_IMPL_PARSE_IPP
#define BOOST_URL_IMPL_PARSE_IPP

#include <boost/url/detail/config.hpp>
#include <boost/url/parse.hpp>
#include <boost/url/rfc/absolute_uri_rule.hpp>
#include <boost/url/rfc/relative_ref_rule.hpp>
#include <boost/url/rfc/uri_rule.hpp>
#include <boost/url/rfc/uri_reference_rule.hpp>
#include <boost/url/rfc/origin_form_rule.hpp>
#include <boost/url/grammar/parse.hpp>

namespace boost {
namespace urls {

system::result<url_view>
parse_absolute_uri(
    core::string_view s)
{
    return grammar::parse(
        s, absolute_uri_rule);
}

system::result<url_view>
parse_origin_form(
    core::string_view s)
{
    return grammar::parse(
        s, origin_form_rule);
}

system::result<url_view>
parse_relative_ref(
    core::string_view s)
{
    return grammar::parse(
        s, relative_ref_rule);
}
system::result<url_view>
parse_uri(
    core::string_view s)
{
    return grammar::parse(
        s, uri_rule);
}

system::result<url_view>
parse_uri_reference(
    core::string_view s)
{
    return grammar::parse(
        s, uri_reference_rule);
}

} // urls
} // boost

#endif
