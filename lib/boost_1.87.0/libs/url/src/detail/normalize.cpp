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
#include <boost/url/decode_view.hpp>
#include "decode.hpp"
#include <boost/url/segments_encoded_view.hpp>
#include <boost/url/grammar/ci_string.hpp>
#include <boost/assert.hpp>
#include <boost/core/ignore_unused.hpp>
#include <cstring>
#include "normalize.hpp"

namespace boost {
namespace urls {
namespace detail {

void
pop_encoded_front(
    core::string_view& s,
    char& c,
    std::size_t& n) noexcept
{
    if(s.front() != '%')
    {
        c = s.front();
        s.remove_prefix(1);
    }
    else
    {
        detail::decode_unsafe(
            &c,
            &c + 1,
            s.substr(0, 3));
        s.remove_prefix(3);
    }
    ++n;
}

int
compare_encoded(
    core::string_view lhs,
    core::string_view rhs) noexcept
{
    std::size_t n0 = 0;
    std::size_t n1 = 0;
    char c0 = 0;
    char c1 = 0;
    while(
        !lhs.empty() &&
        !rhs.empty())
    {
        pop_encoded_front(lhs, c0, n0);
        pop_encoded_front(rhs, c1, n1);
        if (c0 < c1)
            return -1;
        if (c1 < c0)
            return 1;
    }
    n0 += detail::decode_bytes_unsafe(lhs);
    n1 += detail::decode_bytes_unsafe(rhs);
    if (n0 == n1)
        return 0;
    if (n0 < n1)
        return -1;
    return 1;
}

void
digest_encoded(
    core::string_view s,
    fnv_1a& hasher) noexcept
{
    char c = 0;
    std::size_t n = 0;
    while(!s.empty())
    {
        pop_encoded_front(s, c, n);
        hasher.put(c);
    }
}

int
ci_compare_encoded(
    core::string_view lhs,
    core::string_view rhs) noexcept
{
    std::size_t n0 = 0;
    std::size_t n1 = 0;
    char c0 = 0;
    char c1 = 0;
    while (
        !lhs.empty() &&
        !rhs.empty())
    {
        pop_encoded_front(lhs, c0, n0);
        pop_encoded_front(rhs, c1, n1);
        c0 = grammar::to_lower(c0);
        c1 = grammar::to_lower(c1);
        if (c0 < c1)
            return -1;
        if (c1 < c0)
            return 1;
    }
    n0 += detail::decode_bytes_unsafe(lhs);
    n1 += detail::decode_bytes_unsafe(rhs);
    if (n0 == n1)
        return 0;
    if (n0 < n1)
        return -1;
    return 1;
}

void
ci_digest_encoded(
    core::string_view s,
    fnv_1a& hasher) noexcept
{
    char c = 0;
    std::size_t n = 0;
    while(!s.empty())
    {
        pop_encoded_front(s, c, n);
        c = grammar::to_lower(c);
        hasher.put(c);
    }
}

int
compare(
    core::string_view lhs,
    core::string_view rhs) noexcept
{
    auto rlen = (std::min)(lhs.size(), rhs.size());
    for (std::size_t i = 0; i < rlen; ++i)
    {
        char c0 = lhs[i];
        char c1 = rhs[i];
        if (c0 < c1)
            return -1;
        if (c1 < c0)
            return 1;
    }
    if ( lhs.size() == rhs.size() )
        return 0;
    if ( lhs.size() < rhs.size() )
        return -1;
    return 1;
}

int
ci_compare(
    core::string_view lhs,
    core::string_view rhs) noexcept
{
    auto rlen = (std::min)(lhs.size(), rhs.size());
    for (std::size_t i = 0; i < rlen; ++i)
    {
        char c0 = grammar::to_lower(lhs[i]);
        char c1 = grammar::to_lower(rhs[i]);
        if (c0 < c1)
            return -1;
        if (c1 < c0)
            return 1;
    }
    if ( lhs.size() == rhs.size() )
        return 0;
    if ( lhs.size() < rhs.size() )
        return -1;
    return 1;
}

void
ci_digest(
    core::string_view s,
    fnv_1a& hasher) noexcept
{
    for (char c: s)
    {
        c = grammar::to_lower(c);
        hasher.put(c);
    }
}

/* Check if a string ends with the specified suffix (decoded comparison)

   This function determines if a string ends with the specified suffix
   when the string and suffix are compared after percent-decoding.

   @param str The string to check (percent-encoded)
   @param suffix The suffix to check for (percent-decoded)
   @return The number of encoded chars consumed in the string
 */
std::size_t
path_ends_with(
    core::string_view str,
    core::string_view suffix) noexcept
{
    BOOST_ASSERT(!str.empty());
    BOOST_ASSERT(!suffix.empty());
    BOOST_ASSERT(!suffix.contains("%2F"));
    BOOST_ASSERT(!suffix.contains("%2f"));
    auto consume_last = [](
        core::string_view::iterator& it,
        core::string_view::iterator& end,
        char& c)
    {
        BOOST_ASSERT(end > it);
        BOOST_ASSERT(it != end);
        if ((end - it) < 3 ||
            *(std::prev(end, 3)) != '%')
        {
            c = *--end;
            return false;
        }
        detail::decode_unsafe(
            &c,
            &c + 1,
            core::string_view(std::prev(
                end, 3), 3));
        end -= 3;
        return true;
    };

    auto it0 = str.begin();
    auto end0 = str.end();
    auto it1 = suffix.begin();
    auto end1 = suffix.end();
    char c0 = 0;
    char c1 = 0;
    while(
        it0 < end0 &&
        it1 < end1)
    {
        bool const is_encoded = consume_last(it0, end0, c0);
        // The suffix never contains an encoded slash (%2F), and a decoded
        // slash is not equivalent to an encoded slash
        if (is_encoded && c0 == '/')
            return 0;
        consume_last(it1, end1, c1);
        if (c0 != c1)
            return 0;
    }
    bool const consumed_suffix = it1 == end1;
    if (consumed_suffix)
    {
        std::size_t const consumed_encoded = str.end() - end0;
        return consumed_encoded;
    }
    return 0;
}

std::size_t
remove_dot_segments(
    char* dest0,
    char const* end,
    core::string_view input) noexcept
{
    // 1. The input buffer `s` is initialized with
    // the now-appended path components and the
    // output buffer `dest0` is initialized to
    // the empty string.
    char* dest = dest0;
    bool const is_absolute = input.starts_with('/');

    // Step 2 is a loop through 5 production rules:
    // https://www.rfc-editor.org/rfc/rfc3986#section-5.2.4
    //
    // There are no transitions between all rules,
    // which enables some optimizations.
    //
    // Initial:
    // - Rule A: handle initial dots
    // If the input buffer begins with a
    // prefix of "../" or "./", then remove
    // that prefix from the input buffer.
    // Rule A can only happen at the beginning.
    // Errata 4547: Keep "../" in the beginning
    // https://www.rfc-editor.org/errata/eid4547
    //
    // Then:
    // - Rule D: ignore a final ".." or "."
    // if the input buffer consists only  of "."
    // or "..", then remove that from the input
    // buffer.
    // Rule D can only happen after Rule A because:
    // - B and C write "/" to the input
    // - E writes "/" to input or returns
    //
    // Then:
    // - Rule B: ignore ".": write "/" to the input
    // - Rule C: apply "..": remove seg and write "/"
    // - Rule E: copy complete segment
    auto append =
        [](char*& first, char const* last, core::string_view in)
    {
        // append `in` to `dest`
        BOOST_ASSERT(in.size() <= std::size_t(last - first));
        std::memmove(first, in.data(), in.size());
        first += in.size();
        ignore_unused(last);
    };

    auto dot_starts_with = [](
        core::string_view str, core::string_view dots, std::size_t& n)
    {
        // starts_with for encoded/decoded dots
        // or decoded otherwise. return how many
        // chars in str match the dots
        n = 0;
        for (char c: dots)
        {
            if (str.starts_with(c))
            {
                str.remove_prefix(1);
                ++n;
                continue;
            }

            // In the general case, we would need to
            // check if the next char is an encoded
            // dot.
            // However, an encoded dot in `str`
            // would have already been decoded in
            // url_base::normalize_path().
            // This needs to be undone if
            // `remove_dot_segments` is used in a
            // different context.
            // if (str.size() > 2 &&
            //     c == '.'
            //     &&
            //     str[0] == '%' &&
            //     str[1] == '2' &&
            //     (str[2] == 'e' ||
            //      str[2] == 'E'))
            // {
            //     str.remove_prefix(3);
            //     n += 3;
            //     continue;
            // }

            n = 0;
            return false;
        }
        return true;
    };

    auto dot_equal = [&dot_starts_with](
        core::string_view str, core::string_view dots)
    {
        std::size_t n = 0;
        dot_starts_with(str, dots, n);
        return n == str.size();
    };

    // Rule A
    std::size_t n;
    while (!input.empty())
    {
        if (dot_starts_with(input, "../", n))
        {
            // Errata 4547
            append(dest, end, "../");
            input.remove_prefix(n);
            continue;
        }
        else if (!dot_starts_with(input, "./", n))
        {
            break;
        }
        input.remove_prefix(n);
    }

    // Rule D
    if( dot_equal(input, "."))
    {
        input = {};
    }
    else if( dot_equal(input, "..") )
    {
        // Errata 4547
        append(dest, end, "..");
        input = {};
    }

    // 2. While the input buffer is not empty,
    // loop as follows:
    while (!input.empty())
    {
        // Rule B
        bool const is_dot_seg = dot_starts_with(input, "/./", n);
        if (is_dot_seg)
        {
            input.remove_prefix(n - 1);
            continue;
        }

        bool const is_final_dot_seg = dot_equal(input, "/.");
        if (is_final_dot_seg)
        {
            // We can't remove "." from a core::string_view
            // So what we do here is equivalent to
            // replacing s with '/' as required
            // in Rule B and executing the next
            // iteration, which would append this
            // '/' to  the output, as required by
            // Rule E
            append(dest, end, input.substr(0, 1));
            input = {};
            break;
        }

        // Rule C
        bool const is_dotdot_seg = dot_starts_with(input, "/../", n);
        if (is_dotdot_seg)
        {
            core::string_view cur_out(dest0, dest - dest0);
            std::size_t p = cur_out.find_last_of('/');
            bool const has_multiple_segs = p != core::string_view::npos;
            if (has_multiple_segs)
            {
                // output has multiple segments
                // "erase" [p, end] if not "/.."
                core::string_view last_seg(dest0 + p, dest - (dest0 + p));
                bool const prev_is_dotdot_seg = dot_equal(last_seg, "/..");
                if (!prev_is_dotdot_seg)
                {
                    dest = dest0 + p;
                }
                else
                {
                    append(dest, end, "/..");
                }
            }
            else if (dest0 != dest)
            {
                // Only one segment in the output: remove it
                core::string_view last_seg(dest0, dest - dest0);
                bool const prev_is_dotdot_seg = dot_equal(last_seg, "..");
                if (!prev_is_dotdot_seg)
                {
                    dest = dest0;
                    if (!is_absolute)
                    {
                        input.remove_prefix(1);
                    }
                }
                else
                {
                    append(dest, end, "/..");
                }
            }
            else
            {
                // Output is empty
                if (is_absolute)
                {
                    append(dest, end, "/..");
                }
                else
                {
                    // AFREITAS: Although we have no formal proof
                    // for that, the output can't be relative
                    // and empty at this point because relative
                    // paths will fall in the `dest0 != dest`
                    // case above of this rule C and then the
                    // general case of rule E for "..".
                    append(dest, end, "..");
                }
            }
            input.remove_prefix(n - 1);
            continue;
        }

        bool const is_final_dotdot_seg = dot_equal(input, "/..");
        if (is_final_dotdot_seg)
        {
            core::string_view cur_out(dest0, dest - dest0);
            std::size_t p = cur_out.find_last_of('/');
            bool const has_multiple_segs = p != core::string_view::npos;
            if (has_multiple_segs)
            {
                // output has multiple segments
                // "erase" [p, end] if not "/.."
                core::string_view last_seg(dest0 + p, dest - (dest0 + p));
                bool const prev_is_dotdot_seg = dot_equal(last_seg, "/..");
                if (!prev_is_dotdot_seg)
                {
                    dest = dest0 + p;
                    append(dest, end, "/");
                }
                else
                {
                    append(dest, end, "/..");
                }
            }
            else if (dest0 != dest)
            {
                // Only one segment in the output: remove it
                core::string_view last_seg(dest0, dest - dest0);
                bool const prev_is_dotdot_seg = dot_equal(last_seg, "..");
                if (!prev_is_dotdot_seg) {
                    dest = dest0;
                }
                else
                {
                    append(dest, end, "/..");
                }
            }
            else
            {
                // Output is empty: append dotdot
                if (is_absolute)
                {
                    append(dest, end, "/..");
                }
                else
                {
                    // AFREITAS: Although we have no formal proof
                    // for that, the output can't be relative
                    // and empty at this point because relative
                    // paths will fall in the `dest0 != dest`
                    // case above of this rule C and then the
                    // general case of rule E for "..".
                    append(dest, end, "..");
                }
            }
            input = {};
            break;
        }

        // Rule E
        std::size_t p = input.find_first_of('/', 1);
        if (p != core::string_view::npos)
        {
            append(dest, end, input.substr(0, p));
            input.remove_prefix(p);
        }
        else
        {
            append(dest, end, input);
            input = {};
        }
    }

    // 3. Finally, the output buffer is set
    // as the result of remove_dot_segments,
    // and we return its size
    return dest - dest0;
}

char
path_pop_back( core::string_view& s )
{
    if (s.size() < 3 ||
        *std::prev(s.end(), 3) != '%')
    {
        char c = s.back();
        s.remove_suffix(1);
        return c;
    }
    char c = 0;
    detail::decode_unsafe(
        &c, &c + 1, s.substr(s.size() - 3));
    if (c != '/')
    {
        s.remove_suffix(3);
        return c;
    }
    c = s.back();
    s.remove_suffix(1);
    return c;
};

void
pop_last_segment(
    core::string_view& str,
    core::string_view& seg,
    std::size_t& level,
    bool remove_unmatched) noexcept
{
    seg = {};
    std::size_t n = 0;
    while (!str.empty())
    {
        // B.  if the input buffer begins with a
        // prefix of "/./" or "/.", where "." is
        // a complete path segment, then replace
        // that prefix with "/" in the input
        // buffer; otherwise,
        n = detail::path_ends_with(str, "/./");
        if (n)
        {
            seg = str.substr(str.size() - n);
            str.remove_suffix(n);
            continue;
        }
        n = detail::path_ends_with(str, "/.");
        if (n)
        {
            seg = str.substr(str.size() - n, 1);
            str.remove_suffix(n);
            continue;
        }

        // C. if the input buffer begins with a
        // prefix of "/../" or "/..", where ".."
        // is a complete path segment, then
        // replace that prefix with "/" in the
        // input buffer and remove the last
        // segment and its preceding "/"
        // (if any) from the output buffer
        // otherwise,
        n = detail::path_ends_with(str, "/../");
        if (n)
        {
            seg = str.substr(str.size() - n);
            str.remove_suffix(n);
            ++level;
            continue;
        }
        n = detail::path_ends_with(str, "/..");
        if (n)
        {
            seg = str.substr(str.size() - n);
            str.remove_suffix(n);
            ++level;
            continue;
        }

        // E.  move the first path segment in the
        // input buffer to the end of the output
        // buffer, including the initial "/"
        // character (if any) and any subsequent
        // characters up to, but not including,
        // the next "/" character or the end of
        // the input buffer.
        std::size_t p = str.size() > 1
            ? str.find_last_of('/', str.size() - 2)
            : core::string_view::npos;
        if (p != core::string_view::npos)
        {
            seg = str.substr(p + 1);
            str.remove_suffix(seg.size());
        }
        else
        {
            seg = str;
            str = {};
        }

        if (level == 0)
            return;
        if (!str.empty())
            --level;
    }
    // we still need to skip n_skip + 1
    // but the string is empty
    if (remove_unmatched && level)
    {
        seg = "/";
        level = 0;
        return;
    }
    else if (level)
    {
        if (!seg.empty())
        {
            seg = "/../";
        }
        else
        {
            // AFREITAS: this condition
            // is correct, but it might
            // unreachable.
            seg = "/..";
        }
        --level;
        return;
    }
    seg = {};
}

void
normalized_path_digest(
    core::string_view str,
    bool remove_unmatched,
    fnv_1a& hasher) noexcept
{
    core::string_view seg;
    std::size_t level = 0;
    do
    {
        pop_last_segment(
            str, seg, level, remove_unmatched);
        while (!seg.empty())
        {
            char c = path_pop_back(seg);
            hasher.put(c);
        }
    }
    while (!str.empty());
}

// compare segments as if there were a normalized
int
segments_compare(
    segments_encoded_view seg0,
    segments_encoded_view seg1) noexcept
{
    // calculate path size as if it were normalized
    auto normalized_size =
        [](segments_encoded_view seg) -> std::size_t
    {
        if (seg.empty())
            return seg.is_absolute();

        std::size_t n = 0;
        std::size_t skip = 0;
        auto begin = seg.begin();
        auto it = seg.end();
        while (it != begin)
        {
            --it;
            decode_view dseg = **it;
            if (dseg == "..")
                ++skip;
            else if (dseg != ".")
            {
                if (skip)
                    --skip;
                else
                    n += dseg.size() + 1;
            }
        }
        n += skip * 3;
        n -= !seg.is_absolute();
        return n;
    };

    // find the normalized size for the comparison
    std::size_t n0 = normalized_size(seg0);
    std::size_t n1 = normalized_size(seg1);
    std::size_t n00 = n0;
    std::size_t n10 = n1;

    // consume the last char from a segment range
    auto consume_last =
        [](
            std::size_t& n,
            decode_view& dseg,
            segments_encoded_view::iterator& begin,
            segments_encoded_view::iterator& it,
            decode_view::iterator& cit,
            std::size_t& skip,
            bool& at_slash) -> char
    {
        if (cit != dseg.begin())
        {
            // return last char from current segment
            at_slash = false;
            --cit;
            --n;
            return *cit;
        }

        if (!at_slash)
        {
            // current segment dseg is over and
            // previous char was not a slash
            // so we output one
            at_slash = true;
            --n;
            return '/';
        }

        // current segment dseg is over and
        // last char was already the slash
        // between segments, so take the
        // next final segment to consume
        at_slash = false;
        while (cit == dseg.begin())
        {
            // take next segment
            if (it != begin)
                --it;
            else
                break;
            if (**it == "..")
            {
                // skip next if this is ".."
                ++skip;
            }
            else if (**it != ".")
            {
                if (skip)
                {
                    // discount skips
                    --skip;
                }
                else
                {
                    // or update current seg
                    dseg = **it;
                    cit = dseg.end();
                    break;
                }
            }
        }
        // consume from the new current
        // segment
        --n;
        if (cit != dseg.begin())
        {
            // in the general case, we consume
            // one more character from the end
            --cit;
            return *cit;
        }

        // nothing left to consume in the
        // current and new segment
        if (it == begin)
        {
            // if this is the first
            // segment, the segments are
            // over and there can only
            // be repetitions of "../" to
            // output
            return "/.."[n % 3];
        }
        // at other segments, we need
        // a slash to transition to the
        // next segment
        at_slash = true;
        return '/';
    };

    // consume final segments from seg0 that
    // should not influence the comparison
    auto begin0 = seg0.begin();
    auto it0 = seg0.end();
    decode_view dseg0;
    if (it0 != seg0.begin())
    {
        --it0;
        dseg0 = **it0;
    }
    decode_view::iterator cit0 = dseg0.end();
    std::size_t skip0 = 0;
    bool at_slash0 = true;
    while (n0 > n1)
    {
        consume_last(n0, dseg0, begin0, it0, cit0, skip0, at_slash0);
    }

    // consume final segments from seg1 that
    // should not influence the comparison
    auto begin1 = seg1.begin();
    auto it1 = seg1.end();
    decode_view dseg1;
    if (it1 != seg1.begin())
    {
        --it1;
        dseg1 = **it1;
    }
    decode_view::iterator cit1 = dseg1.end();
    std::size_t skip1 = 0;
    bool at_slash1 = true;
    while (n1 > n0)
    {
        consume_last(n1, dseg1, begin1, it1, cit1, skip1, at_slash1);
    }

    int cmp = 0;
    while (n0)
    {
        char c0 = consume_last(
            n0, dseg0, begin0, it0, cit0, skip0, at_slash0);
        char c1 = consume_last(
            n1, dseg1, begin1, it1, cit1, skip1, at_slash1);
        if (c0 < c1)
            cmp = -1;
        else if (c1 < c0)
            cmp = +1;
    }

    if (cmp != 0)
        return cmp;
    if ( n00 == n10 )
        return 0;
    if ( n00 < n10 )
        return -1;
    return 1;
}

} // detail
} // urls
} // boost


