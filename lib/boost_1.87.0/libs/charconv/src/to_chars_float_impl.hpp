// Copyright 2022 Peter Dimov
// Copyright 2023 Matt Borland
// Copyright 2023 Junekey Jeon
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_CHARCONV_DETAIL_TO_CHARS_FLOAT_IMPL_HPP
#define BOOST_CHARCONV_DETAIL_TO_CHARS_FLOAT_IMPL_HPP

#include "float128_impl.hpp"
#include <boost/charconv/detail/apply_sign.hpp>
#include <boost/charconv/detail/integer_search_trees.hpp>
#include <boost/charconv/detail/memcpy.hpp>
#include <boost/charconv/detail/config.hpp>
#include <boost/charconv/detail/dragonbox/floff.hpp>
#include <boost/charconv/detail/bit_layouts.hpp>
#include <boost/charconv/detail/dragonbox/dragonbox.hpp>
#include <boost/charconv/detail/to_chars_integer_impl.hpp>
#include <boost/charconv/detail/to_chars_result.hpp>
#include <boost/charconv/detail/emulated128.hpp>
#include <boost/charconv/detail/fallback_routines.hpp>
#include <boost/charconv/detail/buffer_sizing.hpp>
#include <boost/charconv/config.hpp>
#include <boost/charconv/chars_format.hpp>
#include <system_error>
#include <type_traits>
#include <array>
#include <limits>
#include <utility>
#include <cstring>
#include <cstdio>
#include <cerrno>
#include <cstdint>
#include <climits>
#include <cmath>

#ifdef BOOST_CHARCONV_DEBUG_FIXED
#include <iomanip>
#include <iostream>
#endif

#if (BOOST_CHARCONV_LDBL_BITS == 80 || BOOST_CHARCONV_LDBL_BITS == 128)
#  include <boost/charconv/detail/ryu/ryu_generic_128.hpp>
#  include <boost/charconv/detail/issignaling.hpp>
#endif

namespace boost {
namespace charconv {
namespace detail {

template <typename Real>
inline to_chars_result to_chars_nonfinite(char* first, char* last, Real value, int classification) noexcept;

#if BOOST_CHARCONV_LDBL_BITS == 128 || defined(BOOST_CHARCONV_HAS_STDFLOAT128) || defined(BOOST_CHARCONV_HAS_FLOAT16) || defined(BOOST_CHARCONV_HAS_BRAINFLOAT16)

template <typename Real>
inline to_chars_result to_chars_nonfinite(char* first, char* last, Real value, int classification) noexcept
{
    std::ptrdiff_t offset {};
    std::ptrdiff_t buffer_size = last - first;

    if (classification == FP_NAN)
    {
        bool is_negative = false;
        const bool is_signaling = issignaling(value);

        if (std::signbit(value))
        {
            is_negative = true;
            *first++ = '-';
        }

        if (is_signaling && buffer_size >= (9 + static_cast<int>(is_negative)))
        {
            std::memcpy(first, "nan(snan)", 9);
            offset = 9 + static_cast<int>(is_negative);
        }
        else if (is_negative && buffer_size >= 9)
        {
            std::memcpy(first, "nan(ind)", 8);
            offset = 8;
        }
        else if (!is_negative && !is_signaling && buffer_size >= 3)
        {
            std::memcpy(first, "nan", 3);
            offset = 3;
        }
        else
        {
            // Avoid buffer overflow
            return { first, std::errc::value_too_large };
        }

    }
    else if (classification == FP_INFINITE)
    {
        auto sign_bit_value = std::signbit(value);
        if (sign_bit_value && buffer_size >= 4)
        {
            std::memcpy(first, "-inf", 4);
            offset = 4;
        }
        else if (!sign_bit_value && buffer_size >= 3)
        {
            std::memcpy(first, "inf", 3);
            offset = 3;
        }
        else
        {
            // Avoid buffer overflow
            return { first, std::errc::value_too_large };
        }
    }
    else
    {
        BOOST_UNREACHABLE_RETURN(first);
    }

    return { first + offset, std::errc() };
}

#endif // BOOST_CHARCONV_LDBL_BITS == 128

#ifdef BOOST_CHARCONV_HAS_FLOAT128

// GCC-5 evaluates the following specialization for other types
#if defined(__GNUC__) && __GNUC__ == 5
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wfloat-conversion"
#endif

template <>
inline to_chars_result to_chars_nonfinite<__float128>(char* first, char* last, __float128 value, int classification) noexcept
{
    std::ptrdiff_t offset {};
    std::ptrdiff_t buffer_size = last - first;

    IEEEbinary128 bits;
    std::memcpy(&bits, &value, sizeof(value));
    const bool is_negative = bits.sign;

    if (classification == FP_NAN)
    {
        const bool is_signaling = issignaling(value);

        if (is_negative)
        {
            *first++ = '-';
        }

        if (is_signaling && buffer_size >= (9 + static_cast<int>(is_negative)))
        {
            std::memcpy(first, "nan(snan)", 9);
            offset = 9 + static_cast<int>(is_negative);
        }
        else if (is_negative && buffer_size >= 9)
        {
            std::memcpy(first, "nan(ind)", 8);
            offset = 8;
        }
        else if (!is_negative && !is_signaling && buffer_size >= 3)
        {
            std::memcpy(first, "nan", 3);
            offset = 3;
        }
        else
        {
            // Avoid buffer overflow
            return { first, std::errc::value_too_large };
        }

    }
    else if (classification == FP_INFINITE)
    {
        if (is_negative && buffer_size >= 4)
        {
            std::memcpy(first, "-inf", 4);
            offset = 4;
        }
        else if (!is_negative && buffer_size >= 3)
        {
            std::memcpy(first, "inf", 3);
            offset = 3;
        }
        else
        {
            // Avoid buffer overflow
            return { first, std::errc::value_too_large };
        }
    }
    else
    {
        BOOST_UNREACHABLE_RETURN(first);
    }

    return { first + offset, std::errc() };
}

#if defined(__GNUC__) && __GNUC__ == 5
# pragma GCC diagnostic pop
#endif

#endif // BOOST_CHARCONV_HAS_FLOAT128

template <typename Unsigned_Integer, typename Real, typename std::enable_if<!std::is_same<Unsigned_Integer, uint128>::value, bool>::type = true>
Unsigned_Integer convert_value(Real value) noexcept
{
    Unsigned_Integer temp;
    std::memcpy(&temp, &value, sizeof(Real));
    return temp;
}

template <typename Unsigned_Integer, typename Real, typename std::enable_if<std::is_same<Unsigned_Integer, uint128>::value, bool>::type = true>
Unsigned_Integer convert_value(Real value) noexcept
{
    trivial_uint128 trivial_bits; // NOLINT : Does not need to be init since we memcpy the next step
    std::memcpy(&trivial_bits, &value, sizeof(Real));
    Unsigned_Integer temp {trivial_bits};
    return temp;
}
#ifdef BOOST_MSVC
# pragma warning(push)
# pragma warning(disable: 4127) // Conditional expression is constant (BOOST_IF_CONSTEXPR in pre-C++17 modes)
#endif

template <typename Real, typename Unsigned_Integer>
inline std::uint64_t extract_exp(Real, Unsigned_Integer uint_value, int significand_bits) noexcept
{
    return static_cast<std::uint64_t>(uint_value >> significand_bits);
}

#if BOOST_CHARCONV_LDBL_BITS == 80

// In the 80-bit long double case we need to ignore the pad on top of the high bits
// and the 64-bits of significand in the low bits
// This is easiest to accomplish with a mask rather than all kinds of shifting
template <>
inline std::uint64_t extract_exp<long double, uint128>(long double, uint128 uint_value, int) noexcept
{

    constexpr auto exp_mask = uint128{UINT64_C(0x7FFF), UINT64_C(0)};
    auto exponent = uint_value & exp_mask;
    return exponent.high;
}

#endif

template <typename Real>
to_chars_result to_chars_hex(char* first, char* last, Real value, int precision) noexcept
{
    // Sanity check our bounds
    const std::ptrdiff_t buffer_size = last - first;
    auto real_precision = get_real_precision<Real>(precision);

    if (precision != -1)
    {
        real_precision = precision;
    }

    if (buffer_size < real_precision || first > last)
    {
        return {last, std::errc::value_too_large};
    }

    // Extract the significand and the exponent
    using type_layout =
        #ifdef BOOST_CHARCONV_HAS_FLOAT16
        typename std::conditional<std::is_same<Real, std::float16_t>::value, ieee754_binary16,
        #endif
        #ifdef BOOST_CHARCONV_HAS_BRAINFLOAT16
        typename std::conditional<std::is_same<Real, std::bfloat16_t>::value, brainfloat16,
        #endif

        typename std::conditional<std::is_same<Real, float>::value, ieee754_binary32,
            typename std::conditional<std::is_same<Real, double>::value, ieee754_binary64,
                    #ifdef BOOST_CHARCONV_HAS_FLOAT128
                    typename std::conditional<std::is_same<Real, __float128>::value, ieee754_binary128,
                    #endif
                        #if BOOST_CHARCONV_LDBL_BITS == 128
                        ieee754_binary128
                        #elif BOOST_CHARCONV_LDBL_BITS == 80
                        ieee754_binary80
                        #else
                        ieee754_binary64
                        #endif
                    #ifdef BOOST_CHARCONV_HAS_FLOAT128
                    >::type
                    #endif
            >::type>::type
        #ifdef BOOST_CHARCONV_HAS_FLOAT16
        >::type
        #endif
        #ifdef BOOST_CHARCONV_HAS_BRAINFLOAT16
        >::type
        #endif
        ;

    using Unsigned_Integer =
            typename std::conditional<sizeof(Real) == sizeof(std::uint16_t), std::uint16_t,
            typename std::conditional<std::is_same<Real, float>::value, std::uint32_t,
            typename std::conditional<std::is_same<Real, double>::value, std::uint64_t, uint128>::type>::type>::type;


    Unsigned_Integer uint_value {convert_value<Unsigned_Integer>(value)};

    const Unsigned_Integer denorm_mask = (Unsigned_Integer(1) << (type_layout::significand_bits)) - 1;
    const Unsigned_Integer significand = uint_value & denorm_mask;
    std::uint64_t exponent = extract_exp(value, uint_value, type_layout::significand_bits);

    BOOST_IF_CONSTEXPR (!(std::is_same<Real, float>::value || std::is_same<Real, double>::value
                        #if BOOST_CHARCONV_LDBL_BITS == 80
                        || std::is_same<Real, long double>::value
                        #endif
                        #ifdef BOOST_CHARCONV_HAS_FLOAT16
                        || std::is_same<Real, std::float16_t>::value
                        #endif
                        #ifdef BOOST_CHARCONV_HAS_BRAINFLOAT16
                        || std::is_same<Real, std::bfloat16_t>::value
                        #endif
                        ))
    {
        exponent += 2;
    }

    // Align the significand to the hexit boundaries (i.e. divisible by 4)
    constexpr auto hex_precision =
                                      #ifdef BOOST_CHARCONV_HAS_FLOAT16
                                      std::is_same<Real, std::float16_t>::value ? 3 :
                                      #endif
                                      #ifdef BOOST_CHARCONV_HAS_BRAINFLOAT16
                                      std::is_same<Real, std::bfloat16_t>::value ? 2 :
                                      #endif
                                      std::is_same<Real, float>::value ? 6 :
                                      std::is_same<Real, double>::value ? 13
                                      #if BOOST_CHARCONV_LDBL_BITS == 80
                                      : std::is_same<Real, long double>::value ? 15
                                      #endif
                                      : 28;

    constexpr auto nibble_bits = CHAR_BIT / 2;
    constexpr auto hex_bits = hex_precision * nibble_bits;
    const Unsigned_Integer hex_mask = (static_cast<Unsigned_Integer>(1) << hex_bits) - 1;

    Unsigned_Integer aligned_significand;
    BOOST_IF_CONSTEXPR (std::is_same<Real, float>::value
                        #ifdef BOOST_CHARCONV_HAS_BRAINFLOAT16
                        || std::is_same<Real, std::bfloat16_t>::value
                        #endif
                        )
    {
        aligned_significand = significand << 1;
    }
    #ifdef BOOST_CHARCONV_HAS_FLOAT16
    else BOOST_IF_CONSTEXPR (std::is_same<Real, std::float16_t>::value)
    {
        aligned_significand = significand << 2; // 10/4 = 2.5 so shift two bits
    }
    #endif
    else
    {
        aligned_significand = significand;
    }

    // Adjust the exponent based on the bias as described in IEEE 754
    std::int64_t unbiased_exponent;
    if (exponent == 0 && significand != 0)
    {
        // Subnormal value since we already handled zero
        unbiased_exponent = 1 + type_layout::exponent_bias;
    }
    else
    {
        aligned_significand |= static_cast<Unsigned_Integer>(1) << hex_bits;
        unbiased_exponent = exponent + type_layout::exponent_bias;
    }

    // Bounds check the exponent
    BOOST_IF_CONSTEXPR (std::is_same<Real, float>::value)
    {
        if (unbiased_exponent > 127)
        {
            unbiased_exponent -= 256;
        }
    }
    else BOOST_IF_CONSTEXPR (std::is_same<Real, double>::value)
    {
        if (unbiased_exponent > 1023)
        {
            unbiased_exponent -= 2048;
        }
    }
    else BOOST_IF_CONSTEXPR (std::is_same<Real, long double>::value)
    {
        #if BOOST_CHARCONV_LDBL_BITS == 80
        while (unbiased_exponent > 16383)
        {
            unbiased_exponent -= 32768;
        }
        unbiased_exponent -= 3;
        #else
        while (unbiased_exponent > 16383)
        {
            unbiased_exponent -= 32768;
        }
        #endif
    }
    #ifdef BOOST_CHARCONV_HAS_FLOAT16
    else BOOST_IF_CONSTEXPR (std::is_same<Real, std::float16_t>::value)
    {
        if (unbiased_exponent > 15)
        {
            unbiased_exponent -= 32;
        }
    }
    #endif
    #ifdef BOOST_CHARCONV_HAS_BRAINFLOAT16
    else BOOST_IF_CONSTEXPR (std::is_same<Real, std::bfloat16_t>::value)
    {
        if (unbiased_exponent > 127)
        {
            unbiased_exponent -= 256;
        }
    }
    #endif
    else
    {
        while (unbiased_exponent > 16383)
        {
            unbiased_exponent -= 32768;
        }
    }

    const std::uint32_t abs_unbiased_exponent = unbiased_exponent < 0 ? static_cast<std::uint32_t>(-unbiased_exponent) :
                                                static_cast<std::uint32_t>(unbiased_exponent);

    // Bounds check
    const std::ptrdiff_t total_length = total_buffer_length(real_precision, abs_unbiased_exponent, (value < 0));
    if (total_length > buffer_size)
    {
        return {last, std::errc::value_too_large};
    }

    // Round if required
    if (real_precision < hex_precision)
    {
        const int lost_bits = (hex_precision - real_precision) * nibble_bits;
        const Unsigned_Integer lsb_bit = aligned_significand;
        const Unsigned_Integer round_bit = aligned_significand << 1;
        const Unsigned_Integer tail_bit = round_bit - 1;
        const Unsigned_Integer round = round_bit & (tail_bit | lsb_bit) & (static_cast<Unsigned_Integer>(1) << lost_bits);
        aligned_significand += round;
    }

    // Print the sign
    if (value < 0)
    {
        *first++ = '-';
    }

    // Print the integral part
    #if BOOST_CHARCONV_LDBL_BITS == 80
    std::uint32_t leading_nibble;
    BOOST_IF_CONSTEXPR (std::is_same<Real, long double>::value)
    {
        leading_nibble = static_cast<std::uint32_t>(significand >> hex_bits);
    }
    else
    {
        leading_nibble = static_cast<std::uint32_t>(aligned_significand >> hex_bits);
    }
    #else
    const auto leading_nibble = static_cast<std::uint32_t>(aligned_significand >> hex_bits);
    #endif

    BOOST_CHARCONV_ASSERT(leading_nibble < 16);
    *first++ = digit_table[leading_nibble];
    
    aligned_significand &= hex_mask;

    // Print the fractional part
    if (real_precision > 0)
    {
        *first++ = '.';
        std::int32_t remaining_bits = hex_bits;

        while (true)
        {
            remaining_bits -= nibble_bits;
            const auto current_nibble = static_cast<std::uint32_t>(aligned_significand >> remaining_bits);
            *first++ = digit_table[current_nibble];

            --real_precision;
            if (real_precision == 0)
            {
                break;
            }
            else if (remaining_bits == 0)
            {
                // Do not print trailing zeros with unspecified precision
                if (precision != -1)
                {
                    std::memset(first, '0', static_cast<std::size_t>(real_precision));
                    first += real_precision;
                }
                break;
            }

            // Mask away the hexit we just printed
            aligned_significand &= (static_cast<Unsigned_Integer>(1) << remaining_bits) - 1;
        }
    }

    // Remove any trailing zeros if the precision was unspecified
    if (precision == -1)
    {
        --first;
        while (*first == '0')
        {
            --first;
        }

        // If we have removed everything get rid of the decimal point as well
        if (*first != '.')
        {
            ++first;
        }
    }

    // Print the exponent
    *first++ = 'p';
    if (unbiased_exponent < 0)
    {
        *first++ = '-';
    }
    else
    {
        *first++ = '+';
    }

    return to_chars_int(first, last, abs_unbiased_exponent);
}

#ifdef BOOST_MSVC
# pragma warning(pop)
#endif

template <typename Real>
to_chars_result to_chars_fixed_impl(char* first, char* last, Real value, chars_format fmt = chars_format::general, int precision = -1) noexcept
{
    const std::ptrdiff_t buffer_size = last - first;
    auto real_precision = get_real_precision<Real>(precision);
    if (buffer_size < real_precision || first > last)
    {
        return {last, std::errc::value_too_large};
    }

    auto abs_value = std::abs(value);

    auto value_struct = boost::charconv::detail::to_decimal(value);
    if (value_struct.is_negative)
    {
        *first++ = '-';
    }

    if (value_struct.significand == 0)
    {
        *first++ = '0';
        if (precision > -1)
        {
            *first++ = '.';
            std::memset(first, '0', static_cast<std::size_t>(precision));
            first += precision;
        }

        return {first, std::errc()};
    }

    const int starting_num_digits = num_digits(value_struct.significand);
    int num_dig = starting_num_digits;
    if (precision != -1)
    {
        while (num_dig > precision + 2)
        {
            value_struct.significand /= 10;
            ++value_struct.exponent;
            --num_dig;
        }

        if (num_dig == precision + 2)
        {
            const auto trailing_dig = value_struct.significand % 10;
            value_struct.significand /= 10;
            ++value_struct.exponent;
            --num_dig;

            if (trailing_dig >= 5)
            {
                ++value_struct.significand;
            }
        }

        // In general formatting we remove trailing 0s
        if (fmt == chars_format::general)
        {
            while (value_struct.significand % 10 == 0)
            {
                value_struct.significand /= 10;
                ++value_struct.exponent;
                --num_dig;
            }
        }
    }

    // Make sure the result will fit in the buffer
    const std::ptrdiff_t total_length = total_buffer_length(num_dig, value_struct.exponent, (value < 0));
    if (total_length > buffer_size)
    {
        return {last, std::errc::value_too_large};
    }

    // Insert leading 0s if needed before printing the significand 
    if (abs_value < 1)
    {
        // Additional bounds check for inserted zeros
        if (-value_struct.exponent - starting_num_digits + 2 > (last - first))
        {
            return {last, std::errc::value_too_large};
        }

        std::memcpy(first, "0.", 2U);
        std::memset(first + 2, '0', -value_struct.exponent - starting_num_digits);
        first += 2 - value_struct.exponent - starting_num_digits;
    }

    auto r = to_chars_integer_impl(first, last, value_struct.significand);
    if (r.ec != std::errc())
    {
        return r;
    }

    // Bounds check
    if (abs_value >= 1)
    {
        if (value_struct.exponent < 0 && -value_struct.exponent < buffer_size)
        {
            std::memmove(r.ptr + value_struct.exponent + 1, r.ptr + value_struct.exponent,
                         static_cast<std::size_t>(-value_struct.exponent));
            std::memset(r.ptr + value_struct.exponent, '.', 1);
            ++r.ptr;
        }

        // Add additional zeros as needed
        if (value_struct.exponent > 0)
        {
            const auto zeros_to_append = static_cast<std::size_t>(value_struct.exponent);

            if (zeros_to_append > static_cast<std::size_t>(last - r.ptr))
            {
                return {last, std::errc::value_too_large};
            }

            std::memset(r.ptr, '0', zeros_to_append);
            r.ptr += zeros_to_append;
        }
    }

    return { r.ptr, std::errc() };
}

template <typename Real>
to_chars_result to_chars_float_impl(char* first, char* last, Real value, chars_format fmt, int precision) noexcept
{
    using Unsigned_Integer = typename std::conditional<std::is_same<Real, double>::value, std::uint64_t, std::uint32_t>::type;

    // Sanity check our bounds
    if (first >= last)
    {
        return {last, std::errc::value_too_large};
    }

    auto abs_value = std::abs(value);
    constexpr auto max_fractional_value = std::is_same<Real, double>::value ? static_cast<Real>(1e16) : static_cast<Real>(1e7);
    constexpr auto max_value = static_cast<Real>((std::numeric_limits<Unsigned_Integer>::max)());

    // Unspecified precision so we always go with the shortest representation
    if (precision == -1)
    {
        if (fmt == boost::charconv::chars_format::general)
        {
            if (abs_value >= 1 && abs_value < max_fractional_value)
            {
                return to_chars_fixed_impl(first, last, value, fmt, precision);
            }
            else if (abs_value >= max_fractional_value && abs_value < max_value)
            {
                if (value < 0)
                {
                    *first++ = '-';
                }
                return to_chars_integer_impl(first, last, static_cast<std::uint64_t>(abs_value));
            }
            else
            {
                return boost::charconv::detail::dragonbox_to_chars(value, first, last, fmt);
            }
        }
        else if (fmt == boost::charconv::chars_format::scientific)
        {
            return boost::charconv::detail::dragonbox_to_chars(value, first, last, fmt);
        }
        else if (fmt == boost::charconv::chars_format::fixed)
        {
            return to_chars_fixed_impl(first, last, value, fmt, precision);
        }
    }
    else
    {
        if (fmt != boost::charconv::chars_format::hex)
        {
            if (fmt == boost::charconv::chars_format::general)
            {
                constexpr int max_output_length = std::is_same<Real, double>::value ? 773 : 117;
                constexpr int max_precision = std::is_same<Real, double>::value ? 767 : 112;
                // We remove trailing zeros, so precision > max_precision is same as precision == max_precision.
                if (precision > max_precision)
                {
                    precision = max_precision;
                }
                char temp_buffer[max_output_length];
                auto result = boost::charconv::detail::floff<boost::charconv::detail::main_cache_full,
                                                             boost::charconv::detail::extended_cache_long>(value, precision,
                                                                                                           temp_buffer,
                                                                                                           temp_buffer + max_output_length,
                                                                                                           fmt);
                auto output_size = static_cast<std::size_t>(result.ptr - temp_buffer);
                if (static_cast<std::size_t>(last - first) < output_size)
                {
                    return {last, std::errc::value_too_large};
                }
                std::memcpy(first, temp_buffer, output_size);
                return {first + output_size, std::errc()};
                
            }
            return boost::charconv::detail::floff<boost::charconv::detail::main_cache_full,
                                                  boost::charconv::detail::extended_cache_long>(value, precision,
                                                                                                first, last, fmt);
        }
    }

    const int classification = std::fpclassify(value);
    switch (classification)
    {
        case FP_INFINITE:
        case FP_NAN:
            // The dragonbox impl will return the correct type of NaN
            return boost::charconv::detail::dragonbox_to_chars(value, first, last, chars_format::general);
        case FP_ZERO:
            if (std::signbit(value))
            {
                *first++ = '-';
            }
            std::memcpy(first, "0p+0", 4); // NOLINT : No null terminator is purposeful
            return {first + 4, std::errc()};
        default:
            // Do nothing
            (void)precision;
    }

    // Hex handles both cases already
    return boost::charconv::detail::to_chars_hex(first, last, value, precision);
}

#if (BOOST_CHARCONV_LDBL_BITS == 80 || BOOST_CHARCONV_LDBL_BITS == 128)

template <>
to_chars_result to_chars_float_impl(char* first, char* last, long double value, chars_format fmt, int precision) noexcept
{
    static_assert(std::numeric_limits<long double>::is_iec559, "Long double must be IEEE 754 compliant");

    const auto classification = std::fpclassify(value);
    #if BOOST_CHARCONV_LDBL_BITS == 128
    if (classification == FP_NAN || classification == FP_INFINITE)
    {
        return boost::charconv::detail::to_chars_nonfinite(first, last, value, classification);
    }
    #else
    if (classification == FP_NAN || classification == FP_INFINITE)
    {
        const auto fd128 = boost::charconv::detail::ryu::long_double_to_fd128(value);
        const auto num_chars = boost::charconv::detail::ryu::generic_to_chars(fd128, first, last - first, fmt, precision);

        if (num_chars > 0)
        {
            return { first + num_chars, std::errc() };
        }
        else
        {
            return {last, std::errc::value_too_large};
        }
    }
    #endif

    // Sanity check our bounds
    const std::ptrdiff_t buffer_size = last - first;
    auto real_precision = boost::charconv::detail::get_real_precision<long double>(precision);
    if (buffer_size < real_precision || first > last)
    {
        return {last, std::errc::value_too_large};
    }

    if (fmt == boost::charconv::chars_format::general || fmt == boost::charconv::chars_format::scientific)
    {
        const auto fd128 = boost::charconv::detail::ryu::long_double_to_fd128(value);
        const auto num_chars = boost::charconv::detail::ryu::generic_to_chars(fd128, first, last - first, fmt, precision);

        if (num_chars > 0)
        {
            return { first + num_chars, std::errc() };
        }
    }
    else if (fmt == boost::charconv::chars_format::hex)
    {
        return boost::charconv::detail::to_chars_hex(first, last, value, precision);
    }
    else if (fmt == boost::charconv::chars_format::fixed)
    {
        const auto fd128 = boost::charconv::detail::ryu::long_double_to_fd128(value);
        const auto num_chars = boost::charconv::detail::ryu::generic_to_chars_fixed(fd128, first, last - first, precision);

        if (num_chars > 0)
        {
            return { first + num_chars, std::errc() };
        }
        else if (num_chars == -static_cast<int>(std::errc::value_too_large))
        {
            return { last, std::errc::value_too_large };
        }
    }

    // Fallback to printf methods
    return boost::charconv::detail::to_chars_printf_impl(first, last, value, fmt, precision);
}

#endif

#ifdef BOOST_CHARCONV_HAS_FLOAT128

template <>
to_chars_result to_chars_float_impl(char* first, char* last, __float128 value, chars_format fmt, int precision) noexcept
{
    // Sanity check our bounds
    if (first >= last)
    {
        return {last, std::errc::value_too_large};
    }

    char* const original_first = first;

    if (isnanq(value))
    {
        return boost::charconv::detail::to_chars_nonfinite(first, last, value, FP_NAN);
    }
    else if (isinfq(value))
    {
        return boost::charconv::detail::to_chars_nonfinite(first, last, value, FP_INFINITE);
    }

    // Sanity check our bounds
    const std::ptrdiff_t buffer_size = last - first;
    auto real_precision = boost::charconv::detail::get_real_precision<__float128>(precision);
    if (buffer_size < real_precision || first > last)
    {
        return {last, std::errc::value_too_large};
    }

    if ((fmt == boost::charconv::chars_format::general || fmt == boost::charconv::chars_format::scientific))
    {
        const auto fd128 = boost::charconv::detail::ryu::float128_to_fd128(value);
        const auto num_chars = boost::charconv::detail::ryu::generic_to_chars(fd128, first, last - first, fmt, precision);

        if (num_chars > 0)
        {
            return { first + num_chars, std::errc() };
        }
        else if (num_chars == -1)
        {
            return {last, std::errc::value_too_large};
        }
    }
    else if (fmt == boost::charconv::chars_format::hex)
    {
        return boost::charconv::detail::to_chars_hex(first, last, value, precision);
    }
    else if (fmt == boost::charconv::chars_format::fixed)
    {
        const auto fd128 = boost::charconv::detail::ryu::float128_to_fd128(value);
        const auto num_chars = boost::charconv::detail::ryu::generic_to_chars_fixed(fd128, first, last - first, precision);

        if (num_chars > 0)
        {
            return { first + num_chars, std::errc() };
        }
        else if (num_chars == -static_cast<int>(std::errc::value_too_large))
        {
            return { last, std::errc::value_too_large };
        }
    }

    first = original_first;
    // Fallback to printf
    return boost::charconv::detail::to_chars_printf_impl(first, last, value, fmt, precision);
}

#endif

#if defined(BOOST_CHARCONV_HAS_FLOAT16) || defined(BOOST_CHARCONV_HAS_BRAINFLOAT16)

template <typename T>
to_chars_result to_chars_16_bit_float_impl(char* first, char* last, T value, chars_format fmt, int precision) noexcept
{
    const auto classification = std::fpclassify(value);

    if (classification == FP_NAN || classification == FP_INFINITE)
    {
        return boost::charconv::detail::to_chars_nonfinite(first, last, value, classification);
    }

    // Sanity check our bounds
    const std::ptrdiff_t buffer_size = last - first;
    auto real_precision = boost::charconv::detail::get_real_precision<T>(precision);
    if (buffer_size < real_precision || first > last)
    {
        return {last, std::errc::value_too_large};
    }

    if (fmt == boost::charconv::chars_format::general || fmt == boost::charconv::chars_format::scientific)
    {
        const auto fd128 = boost::charconv::detail::ryu::float16_t_to_fd128(value);
        const auto num_chars = boost::charconv::detail::ryu::generic_to_chars(fd128, first, last - first, fmt, precision);

        if (num_chars > 0)
        {
            return { first + num_chars, std::errc() };
        }
    }
    else if (fmt == boost::charconv::chars_format::hex)
    {
        return boost::charconv::detail::to_chars_hex(first, last, value, precision);
    }
    else if (fmt == boost::charconv::chars_format::fixed)
    {
        const auto fd128 = boost::charconv::detail::ryu::float16_t_to_fd128(value);
        const auto num_chars = boost::charconv::detail::ryu::generic_to_chars_fixed(fd128, first, last - first, precision);

        if (num_chars > 0)
        {
            return { first + num_chars, std::errc() };
        }
        else if (num_chars == -static_cast<int>(std::errc::value_too_large))
        {
            return { last, std::errc::value_too_large };
        }
    }

    // Fallback to printf methods
    return boost::charconv::detail::to_chars_printf_impl(first, last, value, fmt, precision);
}

#endif

} // namespace detail
} // namespace charconv
} // namespace detail

#endif //BOOST_CHARCONV_DETAIL_TO_CHARS_FLOAT_IMPL_HPP
