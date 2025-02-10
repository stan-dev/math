// Copyright 2022 Peter Dimov
// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_CHARCONV_DETAIL_FROM_CHARS_FLOAT_IMPL_HPP
#define BOOST_CHARCONV_DETAIL_FROM_CHARS_FLOAT_IMPL_HPP

#include <boost/charconv/detail/config.hpp>
#include <boost/charconv/detail/from_chars_result.hpp>
#include <boost/charconv/detail/parser.hpp>
#include <boost/charconv/detail/compute_float32.hpp>
#include <boost/charconv/detail/compute_float64.hpp>
#include <boost/charconv/detail/bit_layouts.hpp>
#include <boost/charconv/detail/fallback_routines.hpp>
#include <boost/charconv/chars_format.hpp>
#include <system_error>
#include <cstdlib>
#include <cmath>

namespace boost { namespace charconv { namespace detail {

#ifdef BOOST_MSVC
# pragma warning(push)
# pragma warning(disable: 4244) // Implict converion when BOOST_IF_CONSTEXPR expands to if
#elif defined(__GNUC__) && __GNUC__ >= 5
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wmissing-field-initializers"
# pragma GCC diagnostic ignored "-Wfloat-conversion"
#elif defined(__clang__) && __clang_major__ > 7
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wimplicit-float-conversion"
#elif defined(__clang__) && __clang_major__ <= 7
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wconversion"
#endif

template <typename T>
from_chars_result from_chars_float_impl(const char* first, const char* last, T& value, chars_format fmt) noexcept
{
    bool sign {};
    std::uint64_t significand {};
    std::int64_t  exponent {};

    auto r = boost::charconv::detail::parser(first, last, sign, significand, exponent, fmt);
    if (r.ec == std::errc::value_too_large)
    {
        r.ec = std::errc();
        value = sign ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::infinity();
        return r;
    }
    else if (r.ec == std::errc::not_supported)
    {
        r.ec = std::errc();
        if (significand == 0)
        {
            value = sign ? -std::numeric_limits<T>::quiet_NaN() : std::numeric_limits<long double>::quiet_NaN();
        }
        else
        {
            value = sign ? -std::numeric_limits<T>::signaling_NaN() : std::numeric_limits<T>::signaling_NaN();
        }

        return r;
    }
    if (r.ec != std::errc())
    {
        return r;
    }
    else if (significand == 0)
    {
        value = sign ? static_cast<T>(-0.0L) : static_cast<T>(0.0L);
        return r;
    }
    else if (exponent == -1)
    {
        // A full length significand e.g. -1985444280612224 with a power of -1 sometimes
        // fails in compute_float64 but is trivial to calculate
        // Found investigating GitHub issue #47
        value = (sign ? -static_cast<T>(significand) : static_cast<T>(significand)) / 10;
    }

    bool success {};
    T return_val {};
    BOOST_IF_CONSTEXPR (std::is_same<T, float>::value)
    {
        return_val = compute_float32(exponent, significand, sign, success);
    }
    else
    {
        return_val = compute_float64(exponent, significand, sign, success);
    }

    if (!success)
    {
        if (significand == 1 && exponent == 0)
        {
            value = 1;
            r.ptr = last;
            r.ec = std::errc();
        }
        else
        {
            BOOST_IF_CONSTEXPR (std::is_same<T, float>::value)
            {
                #ifndef __INTEL_LLVM_COMPILER
                if (return_val == HUGE_VALF || return_val == -HUGE_VALF)
                #else
                if (return_val >= std::numeric_limits<T>::max() || return_val <= std::numeric_limits<T>::lowest())
                #endif
                {
                    value = return_val;
                    r.ec = std::errc::result_out_of_range;
                }
                else if (exponent < -46)
                {
                    value = sign ? -0.0F : 0.0;
                    r.ec = std::errc::result_out_of_range;
                }
                else
                {
                    r = from_chars_strtod(first, r.ptr, value);
                }
            }
            else BOOST_IF_CONSTEXPR (std::is_same<T, double>::value)
            {
                #ifndef __INTEL_LLVM_COMPILER
                if (return_val == HUGE_VAL || return_val == -HUGE_VAL)
                #else
                if (return_val >= std::numeric_limits<T>::max() || return_val <= std::numeric_limits<T>::lowest())
                #endif
                {
                    value = return_val;
                    r.ec = std::errc::result_out_of_range;
                }
                else if (exponent < -342)
                {
                    value = sign ? -0.0 : 0.0;
                    r.ec = std::errc::result_out_of_range;
                }
                else
                {
                    r = from_chars_strtod(first, r.ptr, value);
                }
            }
            else BOOST_IF_CONSTEXPR (std::is_same<T, long double>::value)
            {
                #ifndef __INTEL_LLVM_COMPILER
                if (return_val == HUGE_VALL || return_val == -HUGE_VALL)
                #else
                if (return_val >= std::numeric_limits<T>::max() || return_val <= std::numeric_limits<T>::lowest())
                #endif
                {
                    value = return_val;
                    r.ec = std::errc::result_out_of_range;
                }
                #if BOOST_CHARCONV_LDBL_BITS == 64
                else if (exponent < -342)
                #else // 80 or 128 bit long doubles have same range
                else if (exponent < -4965)
                #endif
                {
                    value = sign ? -0.0L : 0.0L;
                    r.ec = std::errc::result_out_of_range;
                }
                else
                {
                    r = from_chars_strtod(first, r.ptr, value);
                }
            }
        }
    }
    else
    {
        value = return_val;
    }

    return r;
}


#ifdef BOOST_MSVC
# pragma warning(pop)
#elif defined(__GNUC__) && __GNUC__ >= 5
# pragma GCC diagnostic pop
#elif defined(__clang__)
# pragma clang diagnostic pop
#endif

}}} // Namespace boost::charconv::detail

#endif // BOOST_CHARCONV_DETAIL_FROM_CHARS_FLOAT_IMPL_HPP
