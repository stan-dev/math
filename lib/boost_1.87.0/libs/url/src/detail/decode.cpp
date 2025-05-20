//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//


#include <boost/url/detail/config.hpp>
#include "decode.hpp"
#include <boost/url/grammar/charset.hpp>
#include <boost/url/grammar/hexdig_chars.hpp>
#include <memory>

namespace boost {
namespace urls {
namespace detail {

char
decode_one(
    char const* const it) noexcept
{
    auto d0 = grammar::hexdig_value(it[0]);
    auto d1 = grammar::hexdig_value(it[1]);
    return static_cast<char>(
        ((static_cast<
            unsigned char>(d0) << 4) +
        (static_cast<
            unsigned char>(d1))));
}

std::size_t
decode_bytes_unsafe(
    core::string_view s) noexcept
{
    auto p = s.begin();
    auto const end = s.end();
    std::size_t dn = 0;
    if(s.size() >= 3)
    {
        auto const safe_end = end - 2;
        while(p < safe_end)
        {
            if(*p != '%')
                p += 1;
            else
                p += 3;
            ++dn;
        }
    }
    dn += end - p;
    return dn;
}

template <bool SpaceAsPlus>
std::size_t
decode_unsafe_is_plus_impl(char c);

template <>
std::size_t
decode_unsafe_is_plus_impl<true>(char c)
{
    return c == '+';
}

template <>
std::size_t
decode_unsafe_is_plus_impl<false>(char)
{
    return false;
}


template <bool SpaceAsPlus>
std::size_t
decode_unsafe_impl(
    char* const dest0,
    char const* end,
    core::string_view s) noexcept
{
    auto it = s.data();
    auto const last = it + s.size();
    auto dest = dest0;

    while(it != last)
    {
        // LCOV_EXCL_START
        if(dest == end)
        {
            /*
             * dest too small: unreachable
             * public functions always pass
             * a buffer of sufficient size
             */
            return dest - dest0;
        }
        // LCOV_EXCL_STOP
        if(decode_unsafe_is_plus_impl<SpaceAsPlus>(*it))
        {
            // plus to space
            *dest++ = ' ';
            ++it;
            continue;
        }
        if(*it == '%')
        {
            // escaped
            ++it;
            // LCOV_EXCL_START
            if(last - it < 2)
            {
                // `%` not followed by two hex digits
                // invalid input: unreachable
                // public functions always pass
                // a valid string_view.
                // initialize output
                std::memset(dest,
                    0, end - dest);
                return dest - dest0;
            }
            // LCOV_EXCL_STOP
            *dest++ = decode_one(it);
            it += 2;
            continue;
        }
        // unescaped
        *dest++ = *it++;
    }
    return dest - dest0;
}

std::size_t
decode_unsafe(
    char* const dest0,
    char const* end,
    core::string_view s,
    encoding_opts opt) noexcept
{
    if(opt.space_as_plus)
    {
        return decode_unsafe_impl<true>(
            dest0, end, s);
    }
    return decode_unsafe_impl<false>(
        dest0, end, s);
}

} // detail
} // urls
} // boost

