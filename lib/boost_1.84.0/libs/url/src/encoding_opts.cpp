//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

#ifndef BOOST_URL_IMPL_ENCONDING_OPTS_IPP
#define BOOST_URL_IMPL_ENCONDING_OPTS_IPP

#include <boost/url/detail/config.hpp>
#include <boost/url/encoding_opts.hpp>

namespace boost {
namespace urls {

encoding_opts::
encoding_opts(
    bool space_as_plus_,
    bool lower_case_,
    bool disallow_null_) noexcept
    : space_as_plus(space_as_plus_)
    , lower_case(lower_case_)
    , disallow_null(disallow_null_)
{}

} // urls
} // boost

#endif
