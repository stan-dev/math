// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_CHARCONV_FLOAT128_IMPL_HPP
#define BOOST_CHARCONV_FLOAT128_IMPL_HPP

#include <boost/charconv/detail/config.hpp>
#include <boost/charconv/detail/ryu/ryu_generic_128.hpp>
#include <boost/charconv/detail/compute_float80.hpp>
#include <boost/charconv/detail/fallback_routines.hpp>
#include <boost/charconv/detail/issignaling.hpp>
#include <boost/charconv/limits.hpp>
#include <system_error>
#include <cstring>
#include <cstdint>

// Only add in float128 support if the build system says it can
#ifdef BOOST_CHARCONV_HAS_QUADMATH

#include <quadmath.h>

#define BOOST_CHARCONV_HAS_FLOAT128

namespace boost {
namespace charconv {

namespace detail {

// --------------------------------------------------------------------------------------------------------------------
// Ryu
// --------------------------------------------------------------------------------------------------------------------


namespace ryu {

inline struct floating_decimal_128 float128_to_fd128(__float128 d) noexcept
{
#ifdef BOOST_CHARCONV_HAS_INT128
    unsigned_128_type bits = 0;
    std::memcpy(&bits, &d, sizeof(__float128));
#else
    trivial_uint128 trivial_bits;
    std::memcpy(&trivial_bits, &d, sizeof(__float128));
    unsigned_128_type bits {trivial_bits};
#endif

    return generic_binary_to_decimal(bits, 112, 15, false);
}

#  ifdef BOOST_CHARCONV_HAS_STDFLOAT128

inline struct floating_decimal_128 stdfloat128_to_fd128(std::float128_t d) noexcept
{
#ifdef BOOST_CHARCONV_HAS_INT128
    unsigned_128_type bits = 0;
    std::memcpy(&bits, &d, sizeof(std::float128_t));
#else
    trivial_uint128 trivial_bits;
    std::memcpy(&trivial_bits, &d, sizeof(std::float128_t));
    unsigned_128_type bits {trivial_bits};
#endif

    return generic_binary_to_decimal(bits, 112, 15, false);
}

#  endif

} // namespace ryu

// --------------------------------------------------------------------------------------------------------------------
// fast_float
// --------------------------------------------------------------------------------------------------------------------

static constexpr __float128 powers_of_tenq[] = {
    1e0Q,  1e1Q,  1e2Q,  1e3Q,  1e4Q,  1e5Q,  1e6Q,
    1e7Q,  1e8Q,  1e9Q,  1e10Q, 1e11Q, 1e12Q, 1e13Q,
    1e14Q, 1e15Q, 1e16Q, 1e17Q, 1e18Q, 1e19Q, 1e20Q,
    1e21Q, 1e22Q, 1e23Q, 1e24Q, 1e25Q, 1e26Q, 1e27Q,
    1e28Q, 1e29Q, 1e30Q, 1e31Q, 1e32Q, 1e33Q, 1e34Q,
    1e35Q, 1e36Q, 1e37Q, 1e38Q, 1e39Q, 1e40Q, 1e41Q,
    1e42Q, 1e43Q, 1e44Q, 1e45Q, 1e46Q, 1e47Q, 1e48Q,
    1e49Q, 1e50Q, 1e51Q, 1e52Q, 1e53Q, 1e54Q, 1e55Q
};

template <typename Unsigned_Integer>
inline __float128 to_float128(Unsigned_Integer w) noexcept
{
    return static_cast<__float128>(w);
}

template <>
inline __float128 to_float128<uint128>(uint128 w) noexcept
{
    return ldexp(static_cast<__float128>(w.high), 64) + static_cast<__float128>(w.low);
}

template <typename Unsigned_Integer, typename ArrayPtr>
inline __float128 fast_path_float128(std::int64_t q, Unsigned_Integer w, bool negative, ArrayPtr table) noexcept
{
    // The general idea is as follows.
    // if 0 <= s <= 2^64 and if 10^0 <= p <= 10^27
    // Both s and p can be represented exactly
    // because of this s*p and s/p will produce
    // correctly rounded values

    auto ld = to_float128(w);

    if (q < 0)
    {
        ld /= table[-q];
    }
    else
    {
        ld *= table[q];
    }

    if (negative)
    {
        ld = -ld;
    }

    return ld;
}

template <typename Unsigned_Integer>
inline __float128 compute_float128(std::int64_t q, Unsigned_Integer w, bool negative, std::errc& success) noexcept
{
    // GLIBC uses 2^-16444 but MPFR uses 2^-16445 as the smallest subnormal value for 80 bit
    // 39 is the max number of digits in an uint128_t
    static constexpr auto smallest_power = -4951 - 39;
    static constexpr auto largest_power = 4932;

    if (-55 <= q && q <= 48 && w <= static_cast<Unsigned_Integer>(1) << 113)
    {
        success = std::errc();
        return fast_path_float128(q, w, negative, powers_of_tenq);
    }

    if (w == 0)
    {
        success = std::errc();
        return negative ? -0.0Q : 0.0Q;
    }
    else if (q > largest_power)
    {
        success = std::errc::result_out_of_range;
        return negative ? -HUGE_VALQ : HUGE_VALQ;
    }
    else if (q < smallest_power)
    {
        success = std::errc::result_out_of_range;
        return negative ? -0.0Q : 0.0Q;
    }

    success = std::errc::not_supported;
    return 0;
}

// --------------------------------------------------------------------------------------------------------------------
// fallback printf
// --------------------------------------------------------------------------------------------------------------------

template <>
inline to_chars_result to_chars_printf_impl<__float128>(char* first, char* last, __float128 value, chars_format fmt, int precision)
{
    // v % + . + num_digits(INT_MAX) + specifier + null terminator
    // 1 + 1 + 10 + 1 + 1
    char format[14] {};
    std::memcpy(format, "%", 1); // NOLINT : No null terminator is purposeful
    std::size_t pos = 1;

    // precision of -1 is unspecified
    if (precision != -1 && fmt != chars_format::fixed)
    {
        format[pos] = '.';
        ++pos;
        const auto unsigned_precision = static_cast<std::uint32_t>(precision);
        if (unsigned_precision < 10)
        {
            boost::charconv::detail::print_1_digit(unsigned_precision, format + pos);
            ++pos;
        }
        else if (unsigned_precision < 100)
        {
            boost::charconv::detail::print_2_digits(unsigned_precision, format + pos);
            pos += 2;
        }
        else
        {
            boost::charconv::detail::to_chars_int(format + pos, format + sizeof(format), precision);
            pos = std::strlen(format);
        }
    }
    else if (fmt == chars_format::fixed)
    {
        // Force 0 decimal places
        std::memcpy(format + pos, ".0", 2); // NOLINT : No null terminator is purposeful
        pos += 2;
    }

    // Add the type identifier
    format[pos] = 'Q';
    ++pos;

    // Add the format character
    switch (fmt)
    {
        case boost::charconv::chars_format::general:
            format[pos] = 'g';
            break;

        case boost::charconv::chars_format::scientific:
            format[pos] = 'e';
            break;

        case boost::charconv::chars_format::fixed:
            format[pos] = 'f';
            break;

        case boost::charconv::chars_format::hex:
            format[pos] = 'a';
            break;
    }

    const auto rv = quadmath_snprintf(first, static_cast<std::size_t>(last - first), format, value);

    if (rv <= 0)
    {
        return {last, static_cast<std::errc>(errno)};
    }

    return {first + rv, std::errc()};
}

// --------------------------------------------------------------------------------------------------------------------
// fallback strtod
// --------------------------------------------------------------------------------------------------------------------

template <>
inline from_chars_result from_chars_strtod_impl<__float128>(const char* first, const char* last, __float128& value, char* buffer) noexcept
{
    // For strto(f/d)
    // Floating point value corresponding to the contents of str on success.
    // If the converted value falls out of range of corresponding return type, range error occurs and HUGE_VAL, HUGE_VALF or HUGE_VALL is returned.
    // If no conversion can be performed, 0 is returned and *str_end is set to str.

    std::memcpy(buffer, first, static_cast<std::size_t>(last - first));
    buffer[last - first] = '\0';
    convert_string_locale(buffer);

    char* str_end;
    __float128 return_value {};
    from_chars_result r {nullptr, std::errc()};

    return_value = strtoflt128(buffer, &str_end);

    if (return_value == HUGE_VALQ)
    {
        r = {last, std::errc::result_out_of_range};
    }

    // Since this is a fallback routine we are safe to check for 0
    if (return_value == 0 && str_end == last)
    {
        r = {first, std::errc::result_out_of_range};
    }

    if (r)
    {
        value = return_value;
        r = {first + (str_end - buffer), std::errc()};
    }

    return r;
}

template <>
inline from_chars_result from_chars_strtod<__float128>(const char* first, const char* last, __float128& value) noexcept
{
    if (last - first < 1024)
    {
        char buffer[1024];
        return from_chars_strtod_impl(first, last, value, buffer);
    }

    // If the string to be parsed does not fit into the 1024 byte static buffer than we have to allocate a buffer.
    // malloc is used here because it does not throw on allocation failure.

    char* buffer = static_cast<char*>(std::malloc(static_cast<std::size_t>(last - first + 1)));
    if (buffer == nullptr)
    {
        return {first, std::errc::not_enough_memory};
    }

    auto r = from_chars_strtod_impl(first, last, value, buffer);
    std::free(buffer);

    return r;
}

// --------------------------------------------------------------------------------------------------------------------
// nans
// --------------------------------------------------------------------------------------------------------------------

struct words
{
#if BOOST_CHARCONV_ENDIAN_LITTLE_BYTE
    std::uint64_t lo;
    std::uint64_t hi;
#else
    std::uint64_t hi;
    std::uint64_t lo;
#endif
};

inline __float128 nans BOOST_PREVENT_MACRO_SUBSTITUTION () noexcept
{
    words bits;
    bits.hi = UINT64_C(0x7FFF400000000000);
    bits.lo = UINT64_C(0);

    __float128 return_val;
    std::memcpy(&return_val, &bits, sizeof(__float128));
    return return_val;
}

inline __float128 nanq BOOST_PREVENT_MACRO_SUBSTITUTION () noexcept
{
    words bits;
    bits.hi = UINT64_C(0x7FFF800000000000);
    bits.lo = UINT64_C(0);

    __float128 return_val;
    std::memcpy(&return_val, &bits, sizeof(__float128));
    return return_val;
}

template <>
inline bool issignaling<__float128> BOOST_PREVENT_MACRO_SUBSTITUTION (__float128 x) noexcept
{
    words bits;
    std::memcpy(&bits, &x, sizeof(__float128));

    std::uint64_t hi_word = bits.hi;
    std::uint64_t lo_word = bits.lo;

    hi_word ^= UINT64_C(0x0000800000000000);
    hi_word |= (lo_word | -lo_word) >> 63;
    return ((hi_word & INT64_MAX) > UINT64_C(0x7FFF800000000000));
}

} //namespace detail
} //namespace charconv
} //namespace boost

#endif //BOOST_CHARCONV_HAS_FLOAT128

#endif //BOOST_CHARCONV_FLOAT128_IMPL_HPP
