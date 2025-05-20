//
// Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//


#include <boost/url/detail/config.hpp>
#include <boost/url/detail/vformat.hpp>
#include "pattern.hpp"

namespace boost {
namespace urls {
namespace detail {

void
vformat_to(
    url_base& u,
    core::string_view fmt,
    detail::format_args args)
{
    parse_pattern(fmt)
        .value().apply(u, args);
}


} // detail
} // urls
} // boost

