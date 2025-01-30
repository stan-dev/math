//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
// Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//


#include <boost/url/detail/config.hpp>
#include <boost/url/segments_view.hpp>
#include <boost/url/parse_path.hpp>

namespace boost {
namespace urls {

segments_view::
segments_view(
    detail::path_ref const& ref) noexcept
    : segments_base(ref)
{
}

segments_view::
segments_view(
    core::string_view s)
    : segments_base(
        parse_path(s).value(
            BOOST_URL_POS))
{
}

} // urls
} // boost

