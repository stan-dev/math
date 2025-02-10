// Copyright 2022 Peter Dimov
// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// https://stackoverflow.com/questions/38060411/visual-studio-2015-wont-suppress-error-c4996
#ifndef _SCL_SECURE_NO_WARNINGS
# define _SCL_SECURE_NO_WARNINGS
#endif
#ifndef NO_WARN_MBCS_MFC_DEPRECATION
# define NO_WARN_MBCS_MFC_DEPRECATION
#endif

#include "float128_impl.hpp"
#include "from_chars_float_impl.hpp"
#include <boost/charconv/detail/fast_float/fast_float.hpp>
#include <boost/charconv/from_chars.hpp>
#include <boost/charconv/detail/bit_layouts.hpp>
#include <system_error>
#include <string>
#include <cstdlib>
#include <cerrno>
#include <cstring>
#include <limits>

#if BOOST_CHARCONV_LDBL_BITS > 64
#  include <boost/charconv/detail/compute_float80.hpp>
#  include <boost/charconv/detail/emulated128.hpp>
#endif

#if defined(__GNUC__) && __GNUC__ < 5
# pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif

boost::charconv::from_chars_result boost::charconv::from_chars_erange(const char* first, const char* last, float& value, boost::charconv::chars_format fmt) noexcept
{
    if (fmt != boost::charconv::chars_format::hex)
    {
        return boost::charconv::detail::fast_float::from_chars(first, last, value, fmt);
    }
    return boost::charconv::detail::from_chars_float_impl(first, last, value, fmt);
}

boost::charconv::from_chars_result boost::charconv::from_chars_erange(const char* first, const char* last, double& value, boost::charconv::chars_format fmt) noexcept
{
    if (fmt != boost::charconv::chars_format::hex)
    {
        return boost::charconv::detail::fast_float::from_chars(first, last, value, fmt);
    }
    return boost::charconv::detail::from_chars_float_impl(first, last, value, fmt);
}

#ifdef BOOST_CHARCONV_HAS_QUADMATH
boost::charconv::from_chars_result boost::charconv::from_chars_erange(const char* first, const char* last, __float128& value, boost::charconv::chars_format fmt) noexcept
{
    bool sign {};
    std::int64_t exponent {};

    #if defined(BOOST_CHARCONV_HAS_INT128) && ((defined(__clang_major__) && __clang_major__ > 12 ) || \
        (defined(BOOST_GCC) && BOOST_GCC > 100000))

    boost::uint128_type significand {};

    #else
    boost::charconv::detail::uint128 significand {};
    #endif

    auto r = boost::charconv::detail::parser(first, last, sign, significand, exponent, fmt);
    if (r.ec == std::errc::value_too_large)
    {
        r.ec = std::errc();

        #if BOOST_CHARCONV_HAS_BUILTIN(__builtin_inf)
        value = sign ? -static_cast<__float128>(__builtin_inf()) : static_cast<__float128>(__builtin_inf());
        #else // Conversion from HUGE_VALL should work
        value = sign ? -static_cast<__float128>(HUGE_VALL) : static_cast<__float128>(HUGE_VALL);
        #endif

        return r;
    }
    else if (r.ec == std::errc::not_supported)
    {
        r.ec = std::errc();
        if (significand == 0)
        {
            #if BOOST_CHARCONV_HAS_BUILTIN(__builtin_nanq)
            value = sign ? -static_cast<__float128>(__builtin_nanq("")) : static_cast<__float128>(__builtin_nanq(""));
            #elif BOOST_CHARCONV_HAS_BUILTIN(__nanq)
            value = sign ? -static_cast<__float128>(__nanq("")) : static_cast<__float128>(__nanq(""));
            #else
            value = boost::charconv::detail::nanq();
            value = sign ? -value : value;
            #endif
        }
        else
        {
            #if BOOST_CHARCONV_HAS_BUILTIN(__builtin_nansq)
            value = sign ? -static_cast<__float128>(__builtin_nansq("")) : static_cast<__float128>(__builtin_nansq(""));
            #elif BOOST_CHARCONV_HAS_BUILTIN(__nansq)
            value = sign ? -static_cast<__float128>(__nansq("")) : static_cast<__float128>(__nansq(""));
            #else
            value = boost::charconv::detail::nans();
            value = sign ? -value : value;
            #endif
        }

        return r;
    }
    else if (r.ec != std::errc())
    {
        return r;
    }
    else if (significand == 0)
    {
        value = sign ? -0.0Q : 0.0Q;
        return r;
    }

    std::errc success {};
    auto return_val = boost::charconv::detail::compute_float128(exponent, significand, sign, success);
    r.ec = static_cast<std::errc>(success);

    if (r.ec == std::errc() || r.ec == std::errc::result_out_of_range)
    {
        value = return_val;
    }
    else if (r.ec == std::errc::not_supported)
    {
        // Fallback routine
        r = boost::charconv::detail::from_chars_strtod(first, last, value);
    }

    return r;
}
#endif

#ifdef BOOST_CHARCONV_HAS_FLOAT16
boost::charconv::from_chars_result boost::charconv::from_chars_erange(const char* first, const char* last, std::float16_t& value, boost::charconv::chars_format fmt) noexcept
{
    float f;
    auto r = boost::charconv::from_chars_erange(first, last, f, fmt);
    if (r.ec == std::errc())
    {
        // Since we are using an interchange format the result could exceed the range of float16_t
        // update the return value or r.ec accordingly
        auto temp = static_cast<std::float16_t>(f);

        if (std::isinf(f) || !std::isinf(temp))
        {
            value = temp;
        }
        else
        {
            r.ec = std::errc::result_out_of_range;
        }
    }

    return r;
}
#endif

#ifdef BOOST_CHARCONV_HAS_FLOAT32
boost::charconv::from_chars_result boost::charconv::from_chars_erange(const char* first, const char* last, std::float32_t& value, boost::charconv::chars_format fmt) noexcept
{
    static_assert(std::numeric_limits<std::float32_t>::digits == FLT_MANT_DIG &&
                  std::numeric_limits<std::float32_t>::min_exponent == FLT_MIN_EXP,
                  "float and std::float32_t are not the same layout like they should be");
    
    float f;
    std::memcpy(&f, &value, sizeof(float));
    const auto r = boost::charconv::from_chars_erange(first, last, f, fmt);
    std::memcpy(&value, &f, sizeof(std::float32_t));
    return r;
}
#endif

#ifdef BOOST_CHARCONV_HAS_FLOAT64
boost::charconv::from_chars_result boost::charconv::from_chars_erange(const char* first, const char* last, std::float64_t& value, boost::charconv::chars_format fmt) noexcept
{
    static_assert(std::numeric_limits<std::float64_t>::digits == DBL_MANT_DIG &&
                  std::numeric_limits<std::float64_t>::min_exponent == DBL_MIN_EXP,
                  "double and std::float64_t are not the same layout like they should be");
    
    double d;
    std::memcpy(&d, &value, sizeof(double));
    const auto r = boost::charconv::from_chars_erange(first, last, d, fmt);
    std::memcpy(&value, &d, sizeof(std::float64_t));
    return r;
}
#endif

#ifdef BOOST_CHARCONV_HAS_BRAINFLOAT16
boost::charconv::from_chars_result boost::charconv::from_chars_erange(const char* first, const char* last, std::bfloat16_t& value, boost::charconv::chars_format fmt) noexcept
{
    float f;
    auto r = boost::charconv::from_chars_erange(first, last, f, fmt);
    if (r.ec == std::errc())
    {
        // Since we are using an interchange format the result could exceed the range of float16_t
        // update the return value or r.ec accordingly
        auto temp = static_cast<std::bfloat16_t>(f);

        if (std::isinf(f) || !std::isinf(temp))
        {
            value = temp;
        }
        else
        {
            r.ec = std::errc::result_out_of_range;
        }
    }

    return r;
}
#endif

#if BOOST_CHARCONV_LDBL_BITS == 64 || defined(BOOST_MSVC)

// Since long double is just a double we use the double implementation and cast into value
boost::charconv::from_chars_result boost::charconv::from_chars_erange(const char* first, const char* last, long double& value, boost::charconv::chars_format fmt) noexcept
{
    static_assert(sizeof(double) == sizeof(long double), "64 bit long double detected, but the size is incorrect");
    
    double d;
    std::memcpy(&d, &value, sizeof(double));
    const auto r = boost::charconv::from_chars_erange(first, last, d, fmt);
    std::memcpy(&value, &d, sizeof(long double));

    return r;
}

#elif !defined(BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE)

boost::charconv::from_chars_result boost::charconv::from_chars_erange(const char* first, const char* last, long double& value, boost::charconv::chars_format fmt) noexcept
{
    static_assert(std::numeric_limits<long double>::is_iec559, "Long double must be IEEE 754 compliant");

    bool sign {};
    std::int64_t exponent {};

    #if defined(BOOST_CHARCONV_HAS_INT128) && ((defined(__clang_major__) && __clang_major__ > 12 ) || \
        (defined(BOOST_GCC) && BOOST_GCC > 100000))

    boost::uint128_type significand {};

    #else
    boost::charconv::detail::uint128 significand {};
    #endif

    auto r = boost::charconv::detail::parser(first, last, sign, significand, exponent, fmt);
    if (r.ec == std::errc::value_too_large)
    {
        r.ec = std::errc();
        value = sign ? -std::numeric_limits<long double>::infinity() : std::numeric_limits<long double>::infinity();
        return r;
    }
    else if (r.ec == std::errc::not_supported)
    {
        r.ec = std::errc();
        if (significand == 0)
        {
            value = sign ? -std::numeric_limits<long double>::quiet_NaN() : std::numeric_limits<long double>::quiet_NaN();
        }
        else
        {
            value = sign ? -std::numeric_limits<long double>::signaling_NaN() : std::numeric_limits<long double>::signaling_NaN();
        }

        return r;
    }
    else if (r.ec != std::errc())
    {
        return r;
    }
    else if (significand == 0)
    {
        value = sign ? -0.0L : 0.0L;
        return r;
    }

    std::errc success {};
    auto return_val = boost::charconv::detail::compute_float80<long double>(exponent, significand, sign, success);
    r.ec = success;

    if (r.ec == std::errc() || r.ec == std::errc::result_out_of_range)
    {
        value = return_val;
    }
    else if (r.ec == std::errc::not_supported)
    {
        // Fallback routine
        r = boost::charconv::detail::from_chars_strtod(first, last, value);
    }

    return r;
}

#if defined(BOOST_CHARCONV_HAS_STDFLOAT128) && defined(BOOST_CHARCONV_HAS_QUADMATH)
boost::charconv::from_chars_result boost::charconv::from_chars_erange(const char* first, const char* last, std::float128_t& value, boost::charconv::chars_format fmt) noexcept
{
    static_assert(sizeof(__float128) == sizeof(std::float128_t));

    __float128 q;
    std::memcpy(&q, &value, sizeof(__float128));
    const auto r = boost::charconv::from_chars_erange(first, last, q, fmt);
    std::memcpy(&value, &q, sizeof(std::float128_t));

    return r;
}
#endif

#endif // long double implementations

// String view overloads

boost::charconv::from_chars_result boost::charconv::from_chars_erange(boost::core::string_view sv, float& value, boost::charconv::chars_format fmt) noexcept
{
    return boost::charconv::from_chars_erange(sv.data(), sv.data() + sv.size(), value, fmt);
}

boost::charconv::from_chars_result boost::charconv::from_chars_erange(boost::core::string_view sv, double & value, boost::charconv::chars_format fmt) noexcept
{
    return boost::charconv::from_chars_erange(sv.data(), sv.data() + sv.size(), value, fmt);
}

#ifndef BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE
boost::charconv::from_chars_result boost::charconv::from_chars_erange(boost::core::string_view sv, long double& value, boost::charconv::chars_format fmt) noexcept
{
    return boost::charconv::from_chars_erange(sv.data(), sv.data() + sv.size(), value, fmt);
}
#endif

#ifdef BOOST_CHARCONV_HAS_QUADMATH
boost::charconv::from_chars_result boost::charconv::from_chars_erange(boost::core::string_view sv, __float128& value, boost::charconv::chars_format fmt) noexcept
{
    return boost::charconv::from_chars_erange(sv.data(), sv.data() + sv.size(), value, fmt);
}
#endif

// <stdfloat> types
#ifdef BOOST_CHARCONV_HAS_FLOAT16
boost::charconv::from_chars_result boost::charconv::from_chars_erange(boost::core::string_view sv, std::float16_t& value, boost::charconv::chars_format fmt) noexcept
{
    return boost::charconv::from_chars_erange(sv.data(), sv.data() + sv.size(), value, fmt);
}
#endif
#ifdef BOOST_CHARCONV_HAS_FLOAT32
boost::charconv::from_chars_result boost::charconv::from_chars_erange(boost::core::string_view sv, std::float32_t& value, boost::charconv::chars_format fmt) noexcept
{
    return boost::charconv::from_chars_erange(sv.data(), sv.data() + sv.size(), value, fmt);
}
#endif
#ifdef BOOST_CHARCONV_HAS_FLOAT64
boost::charconv::from_chars_result boost::charconv::from_chars_erange(boost::core::string_view sv, std::float64_t& value, boost::charconv::chars_format fmt) noexcept
{
    return boost::charconv::from_chars_erange(sv.data(), sv.data() + sv.size(), value, fmt);
}
#endif
#if defined(BOOST_CHARCONV_HAS_STDFLOAT128) && defined(BOOST_CHARCONV_HAS_QUADMATH)
boost::charconv::from_chars_result boost::charconv::from_chars_erange(boost::core::string_view sv, std::float128_t& value, boost::charconv::chars_format fmt) noexcept
{
    return boost::charconv::from_chars_erange(sv.data(), sv.data() + sv.size(), value, fmt);
}
#endif
#ifdef BOOST_CHARCONV_HAS_BRAINFLOAT16
boost::charconv::from_chars_result boost::charconv::from_chars_erange(boost::core::string_view sv, std::bfloat16_t& value, boost::charconv::chars_format fmt) noexcept
{
    return boost::charconv::from_chars_erange(sv.data(), sv.data() + sv.size(), value, fmt);
}
#endif

namespace {

// Adheres to the STL strictly as opposed to fixing the ERANGE problem (which pre-review was the library default behavior)
template <typename T>
boost::charconv::from_chars_result from_chars_strict_impl(const char* first, const char* last, T& value, boost::charconv::chars_format fmt) noexcept
{
    T temp_value {};
    const auto r = boost::charconv::from_chars_erange(first, last, temp_value, fmt);

    if (r)
    {
        value = temp_value;
    }

    return r;
}

}

boost::charconv::from_chars_result boost::charconv::from_chars(const char* first, const char* last, float& value, boost::charconv::chars_format fmt) noexcept
{
    return from_chars_strict_impl(first, last, value, fmt);
}

boost::charconv::from_chars_result boost::charconv::from_chars(const char* first, const char* last, double& value, boost::charconv::chars_format fmt) noexcept
{
    return from_chars_strict_impl(first, last, value, fmt);
}

#ifndef BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE
boost::charconv::from_chars_result boost::charconv::from_chars(const char* first, const char* last, long double& value, boost::charconv::chars_format fmt) noexcept
{
    return from_chars_strict_impl(first, last, value, fmt);
}
#endif

#ifdef BOOST_CHARCONV_HAS_QUADMATH
boost::charconv::from_chars_result boost::charconv::from_chars(const char* first, const char* last, __float128& value, boost::charconv::chars_format fmt) noexcept
{
    return from_chars_strict_impl(first, last, value, fmt);
}
#endif

#ifdef BOOST_CHARCONV_HAS_FLOAT16
boost::charconv::from_chars_result boost::charconv::from_chars(const char* first, const char* last, std::float16_t& value, boost::charconv::chars_format fmt) noexcept
{
    return from_chars_strict_impl(first, last, value, fmt);
}
#endif

#ifdef BOOST_CHARCONV_HAS_FLOAT32
boost::charconv::from_chars_result boost::charconv::from_chars(const char* first, const char* last, std::float32_t& value, boost::charconv::chars_format fmt) noexcept
{
    return from_chars_strict_impl(first, last, value, fmt);
}
#endif

#ifdef BOOST_CHARCONV_HAS_FLOAT64
boost::charconv::from_chars_result boost::charconv::from_chars(const char* first, const char* last, std::float64_t& value, boost::charconv::chars_format fmt) noexcept
{
    return from_chars_strict_impl(first, last, value, fmt);
}
#endif

#if defined(BOOST_CHARCONV_HAS_STDFLOAT128) && defined(BOOST_CHARCONV_HAS_QUADMATH)
boost::charconv::from_chars_result boost::charconv::from_chars(const char* first, const char* last, std::float128_t& value, boost::charconv::chars_format fmt) noexcept
{
    return from_chars_strict_impl(first, last, value, fmt);
}
#endif

#ifdef BOOST_CHARCONV_HAS_BRAINFLOAT16
boost::charconv::from_chars_result boost::charconv::from_chars(const char* first, const char* last, std::bfloat16_t& value, boost::charconv::chars_format fmt) noexcept
{
    return from_chars_strict_impl(first, last, value, fmt);
}
#endif

boost::charconv::from_chars_result boost::charconv::from_chars(boost::core::string_view sv, float& value, boost::charconv::chars_format fmt) noexcept
{
    return from_chars_strict_impl(sv.data(), sv.data() + sv.size(), value, fmt);
}

boost::charconv::from_chars_result boost::charconv::from_chars(boost::core::string_view sv, double& value, boost::charconv::chars_format fmt) noexcept
{
    return from_chars_strict_impl(sv.data(), sv.data() + sv.size(), value, fmt);
}

#ifndef BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE
boost::charconv::from_chars_result boost::charconv::from_chars(boost::core::string_view sv, long double& value, boost::charconv::chars_format fmt) noexcept
{
    return from_chars_strict_impl(sv.data(), sv.data() + sv.size(), value, fmt);
}
#endif

#ifdef BOOST_CHARCONV_HAS_QUADMATH
boost::charconv::from_chars_result boost::charconv::from_chars(boost::core::string_view sv, __float128& value, boost::charconv::chars_format fmt) noexcept
{
    return from_chars_strict_impl(sv.data(), sv.data() + sv.size(), value, fmt);
}
#endif

#ifdef BOOST_CHARCONV_HAS_FLOAT16
boost::charconv::from_chars_result boost::charconv::from_chars(boost::core::string_view sv, std::float16_t& value, boost::charconv::chars_format fmt) noexcept
{
    return from_chars_strict_impl(sv.data(), sv.data() + sv.size(), value, fmt);
}
#endif

#ifdef BOOST_CHARCONV_HAS_FLOAT32
boost::charconv::from_chars_result boost::charconv::from_chars(boost::core::string_view sv, std::float32_t& value, boost::charconv::chars_format fmt) noexcept
{
    return from_chars_strict_impl(sv.data(), sv.data() + sv.size(), value, fmt);
}
#endif

#ifdef BOOST_CHARCONV_HAS_FLOAT64
boost::charconv::from_chars_result boost::charconv::from_chars(boost::core::string_view sv, std::float64_t& value, boost::charconv::chars_format fmt) noexcept
{
    return from_chars_strict_impl(sv.data(), sv.data() + sv.size(), value, fmt);
}
#endif

#if defined(BOOST_CHARCONV_HAS_STDFLOAT128) && defined(BOOST_CHARCONV_HAS_QUADMATH)
boost::charconv::from_chars_result boost::charconv::from_chars(boost::core::string_view sv, std::float128_t& value, boost::charconv::chars_format fmt) noexcept
{
    return from_chars_strict_impl(sv.data(), sv.data() + sv.size(), value, fmt);
}
#endif

#ifdef BOOST_CHARCONV_HAS_BRAINFLOAT16
boost::charconv::from_chars_result boost::charconv::from_chars(boost::core::string_view sv, std::bfloat16_t& value, boost::charconv::chars_format fmt) noexcept
{
    return from_chars_strict_impl(sv.data(), sv.data() + sv.size(), value, fmt);
}
#endif
