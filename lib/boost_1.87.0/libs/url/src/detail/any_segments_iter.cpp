//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//


#include <boost/url/detail/config.hpp>
#include "../rfc/detail/charsets.hpp"
#include <boost/url/detail/any_segments_iter.hpp>
#include <boost/core/detail/string_view.hpp>
#include <boost/url/encode.hpp>
#include <boost/url/rfc/pchars.hpp>

namespace boost {
namespace urls {
namespace detail {

//------------------------------------------------
//
// segment_iter
//
//------------------------------------------------

segment_iter::
segment_iter(
    core::string_view s_) noexcept
    : any_segments_iter(s_)
{
    front = s;
    fast_nseg = 1;
}

void
segment_iter::
rewind() noexcept
{
    at_end_ = false;
}

bool
segment_iter::
measure(
    std::size_t& n) noexcept
{
    if(at_end_)
        return false;
    encoding_opts opt;
    opt.space_as_plus = false;
    n += encoded_size(
        s,
        encode_colons ?
            nocolon_pchars :
            pchars,
        opt);
    at_end_ = true;
    return true;
}

void
segment_iter::
copy(
    char*& dest,
    char const* end) noexcept
{
    encoding_opts opt;
    opt.space_as_plus = false;
    dest += encode(
        dest,
        end - dest,
        s,
        encode_colons ?
            nocolon_pchars :
            pchars,
        opt);
}

//------------------------------------------------
//
// segments_iter_base
//
//------------------------------------------------

void
segments_iter_base::
measure_impl(
    std::size_t& n,
    core::string_view s,
    bool encode_colons) noexcept
{
    encoding_opts opt;
    opt.space_as_plus = false;
    n += encoded_size(
        s,
        encode_colons ?
            nocolon_pchars :
            pchars,
        opt);
}

void
segments_iter_base::
copy_impl(
    char*& dest,
    char const* end,
    core::string_view s,
    bool encode_colons) noexcept
{
    encoding_opts opt;
    opt.space_as_plus = false;
    dest += encode(
        dest,
        end - dest,
        s,
        encode_colons ?
            nocolon_pchars :
            pchars,
        opt);
}

//------------------------------------------------
//
// segment_encoded_iter
//
//------------------------------------------------

segment_encoded_iter::
segment_encoded_iter(
    pct_string_view const& s_) noexcept
    : any_segments_iter(s_)
{
    front = s;
    fast_nseg = 1;
}

void
segment_encoded_iter::
rewind() noexcept
{
    at_end_ = false;
}

bool
segment_encoded_iter::
measure(
    std::size_t& n) noexcept
{
    if(at_end_)
        return false;
    n += detail::re_encoded_size_unsafe(
        s,
        encode_colons ?
            nocolon_pchars :
            pchars);
    at_end_ = true;
    return true;
}

void
segment_encoded_iter::
copy(
    char*& dest,
    char const* end) noexcept
{
    detail::re_encode_unsafe(
        dest,
        end,
        s,
        encode_colons ?
            nocolon_pchars :
            pchars);
}

//------------------------------------------------
//
// segments_encoded_iter_base
//
//------------------------------------------------

void
segments_encoded_iter_base::
measure_impl(
    std::size_t& n,
    core::string_view s,
    bool encode_colons) noexcept
{
    n += detail::re_encoded_size_unsafe(
        s,
        encode_colons ?
            nocolon_pchars :
            pchars);
}

void
segments_encoded_iter_base::
copy_impl(
    char*& dest,
    char const* end,
    core::string_view s,
    bool encode_colons) noexcept
{
    detail::re_encode_unsafe(
        dest,
        end,
        s,
        encode_colons ?
            nocolon_pchars :
            pchars);
}

} // detail
} // urls
} // boost

