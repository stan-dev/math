//  (C) Copyright Stephan T. Lavavej 2019 - 2023.
//  (C) Copyright Matt Borland 2023.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Derived from mailing list post: https://lists.boost.org/Archives/boost/2023/05/254660.php

#include <iostream>

#ifdef BOOST_CHARCONV_RUN_BENCHMARKS
#ifdef __GNUC__
#  pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

#ifndef _MSC_VER
#define AVOID_SPRINTF_S
#endif // _MSC_VER

#ifndef AVOID_CHARCONV
#include <charconv>
#endif // AVOID_CHARCONV
#include <array>
#include <chrono>
#include <exception>
#include <random>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <system_error>
#include <type_traits>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <cmath>

#include <fmt/format.h>
#include <double-conversion/double-conversion.h>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/phoenix/stl.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/charconv.hpp>

using namespace std;
using namespace std::chrono;

void verify(const bool b) {
   if (!b) {
       puts("VERIFICATION FAILURE");
       exit(EXIT_FAILURE);
   }
}

template <typename T>
void verify_fp(T x, T y, int tol = 10) {
    if (abs(boost::math::float_distance(x, y)) > tol)
    {
        puts("VERIFICATION FAILURE");
        exit(EXIT_FAILURE);
        // cerr << "Float distance between: " << x << " and " << y << " is " << abs(boost::math::float_distance(x, y)) << endl;
    }
}

enum class RoundTrip { Sci, Fix, Gen, Hex, Lossy, u32, u64 };

constexpr size_t N = 2000; // how many values to test

constexpr size_t K = 5; // how many times to repeat the test, for cleaner timing

constexpr size_t BufSize = 2'000'000; // more than enough

unsigned int global_dummy = 0;

template <typename Num>
struct scientific_policy : boost::spirit::karma::real_policies<Num>
{
    static bool trailing_zeros(Num)
    {
        return false;
    }

    static unsigned precision(Num)
    {
        return std::numeric_limits<Num>::digits10 + 1;
    }
};

template <typename T>
void test_boost_spirit_karma(const char* const str, const vector<T>& vec) {
    namespace karma = boost::spirit::karma;

    using science_type_double = karma::real_generator<double, scientific_policy<double>>;
    science_type_double const scientific_double = science_type_double();
    using science_type_float = karma::real_generator<float, scientific_policy<float>>;
    science_type_float const scientific_float = science_type_float();

    const auto start = steady_clock::now();

    for (size_t k = 0; k < K; ++k) {
        for (const auto& elem : vec) {
            char buffer[BufSize];
            char* it = buffer;
            
            if constexpr (std::is_same_v<T, float>) {
                karma::generate(it, scientific_float, elem);
            } else if constexpr (std::is_same_v<T, double>) {
                karma::generate(it, scientific_double, elem);
            } else if constexpr (std::is_same_v<T, uint32_t>) {
                karma::generate(it, karma::ulong_, (unsigned long)elem);
            } else if constexpr (std::is_same_v<T, uint64_t>) {
                karma::generate(it, karma::ulong_long, (unsigned long long)elem);
            }
        }
    }
    const auto finish = steady_clock::now();

    printf("%6.1f ns | %s\n", duration<double, nano>{finish - start}.count() / (N * K), str);

    for (const auto& elem : vec) {
        char buffer[BufSize] {};
        char* it = buffer;

        if constexpr (std::is_same_v<T, float>) {
            karma::generate(it, scientific_float, elem);
        } else if constexpr (std::is_same_v<T, double>) {
            karma::generate(it, scientific_double, elem);
        } else if constexpr (std::is_same_v<T, uint32_t>) {
            karma::generate(it, karma::ulong_, (unsigned long)elem);
        } else if constexpr (std::is_same_v<T, uint64_t>) {
            karma::generate(it, karma::ulong_long, (unsigned long long)elem);
        }

        T round_trip;
        auto r = from_chars(buffer, buffer + BufSize, round_trip);
        verify(r.ec == std::errc());

        if constexpr (std::is_integral_v<T>) {
            verify(round_trip == elem);
        } else {
            // Have to set a loose tolerance,
            // because even with the user specified policy rounding errors can be large
            //
            // Float distance between: -0.0036966 and -0.00369657 is 149
            // Float distance between: 0.002663 and 0.00266298 is 88
            // Float distance between: -0.0025576 and -0.00255758 is 107
            verify_fp(round_trip, elem, 500);
        }

    }
}

template <typename T>
void test_lexical_cast(const char* const str, const vector<T>& vec) {
    const auto start = steady_clock::now();
    for (size_t k = 0; k < K; ++k) {
        for (const auto& elem : vec) {
            const auto ret = boost::lexical_cast<array<char, 25>>(elem); // "-1.2345678901234567e-100" plus null term
            global_dummy += static_cast<unsigned int>(ret[0]);
        }
    }
    const auto finish = steady_clock::now();

    printf("%6.1f ns | %s\n", duration<double, nano>{finish - start}.count() / (N * K), str);

    for (const auto& elem : vec) {
        const auto ret = boost::lexical_cast<array<char, 25>>(elem);

        T round_trip;
        auto r = from_chars(ret.data(), ret.data() + ret.size(), round_trip);
        verify(r.ec == std::errc());

        if constexpr (std::is_integral_v<T>) {
            verify(round_trip == elem);
        } else {
            verify_fp(round_trip, elem, 2);
        }

    }
}

template <typename Floating>
int sprintf_wrapper(char (&buf)[BufSize], const char* const fmt, const Floating elem) {
#ifdef AVOID_SPRINTF_S
    return snprintf(buf, sizeof(buf), fmt, elem);
#else // AVOID_SPRINTF_S
    return sprintf_s(buf, BufSize, fmt, elem);
#endif // AVOID_SPRINTF_S
}

template <RoundTrip RT, typename T>
void test_sprintf(const char* const str, const vector<T>& vec, const char* const fmt) {

    char buf[BufSize];

    const auto start = steady_clock::now();
    for (size_t k = 0; k < K; ++k) {
        for (const auto& elem : vec) {
            int ret;
            if constexpr (std::is_same_v<T, uint32_t>)
                ret = sprintf_wrapper(buf, fmt, (unsigned long)elem);
            else if constexpr (std::is_same_v<T, uint64_t>)
                ret = sprintf_wrapper(buf, fmt, (unsigned long long)elem);
            else
                ret = sprintf_wrapper(buf, fmt, elem);

            global_dummy += static_cast<unsigned int>(ret);
            global_dummy += static_cast<unsigned int>(buf[0]);
        }
    }
    const auto finish = steady_clock::now();

    printf("%6.1f ns | %s\n", duration<double, nano>{finish - start}.count() / (N * K), str);

    for (const auto& elem : vec) {
        if constexpr (std::is_same_v<T, uint32_t>)
            verify(sprintf_wrapper(buf, fmt, (unsigned long)elem) != -1);
        else if constexpr (std::is_same_v<T, uint64_t>)
            verify(sprintf_wrapper(buf, fmt, (unsigned long long)elem) != -1);
        else
            verify(sprintf_wrapper(buf, fmt, elem) != -1);

        if constexpr (RT == RoundTrip::Lossy) {
            // skip lossy conversions
        } else if constexpr (is_same_v<T, float>) {
            verify(strtof(buf, nullptr) == elem);
        } else if constexpr (is_same_v<T, double>) {
            verify(strtod(buf, nullptr) == elem);
        } else if constexpr (is_same_v<T, uint32_t>) {
            verify(strtoul(buf, nullptr, 10) == (unsigned long)elem);
        } else if constexpr (is_same_v<T, uint64_t>) {
            verify(strtoull(buf, nullptr, 10) == (unsigned long long)elem);
        }
    }
}

#ifndef AVOID_CHARCONV
constexpr chars_format chars_format_from_RoundTrip(const RoundTrip rt) {
    switch (rt) {
    case RoundTrip::Sci:
        return chars_format::scientific;
    case RoundTrip::Fix:
        return chars_format::fixed;
    case RoundTrip::Gen:
        return chars_format::general;
    case RoundTrip::Hex:
        return chars_format::hex;
    case RoundTrip::Lossy:
    default:
        puts("CHARS FORMAT FAIL");
        exit(EXIT_FAILURE);
    }
}

constexpr boost::charconv::chars_format boost_chars_format_from_RoundTrip(const RoundTrip rt) {
    switch (rt) {
    case RoundTrip::Sci:
        return boost::charconv::chars_format::scientific;
    case RoundTrip::Fix:
        return boost::charconv::chars_format::fixed;
    case RoundTrip::Gen:
        return boost::charconv::chars_format::general;
    case RoundTrip::Hex:
        return boost::charconv::chars_format::hex;
    case RoundTrip::Lossy:
    default:
        puts("CHARS FORMAT FAIL");
        exit(EXIT_FAILURE);
    }
}

template <RoundTrip RT, typename T, typename... Args>
void test_STL_to_chars(const char* const str, const vector<T>& vec, const Args&... args) {

    char buf[BufSize];

    const auto start = steady_clock::now();
    for (size_t k = 0; k < K; ++k) {
        for (const auto& elem : vec) {
            const auto result = to_chars(buf, buf + BufSize, elem, args...);

            global_dummy += static_cast<unsigned int>(result.ptr - buf);
            global_dummy += static_cast<unsigned int>(buf[0]);
        }
    }
    const auto finish = steady_clock::now();

    printf("%6.1f ns | %s\n", duration<double, nano>{finish - start}.count() / (N * K), str);

    for (const auto& elem : vec) {
        const auto result = to_chars(buf, buf + BufSize, elem, args...);
        verify(result.ec == errc{});

        if constexpr (RT == RoundTrip::Lossy) {
            // skip lossy conversions
        } else {
            T round_trip;
            from_chars_result from_result;
            if constexpr (std::is_same_v<T, uint32_t> || std::is_same_v<T, uint64_t>)
                from_result = std::from_chars(buf, result.ptr, round_trip, args...);
            else
                from_result = std::from_chars(buf, result.ptr, round_trip, chars_format_from_RoundTrip(RT));
            verify(from_result.ec == errc{});
            verify(from_result.ptr == result.ptr);
            verify(round_trip == elem);
        }
    }
}

template <RoundTrip RT, typename T, typename... Args>
void test_boost_to_chars(const char* const str, const vector<T>& vec, const Args&... args) {

    char buf[BufSize] {};

    const auto start = steady_clock::now();
    for (size_t k = 0; k < K; ++k) {
        for (const auto& elem : vec) {
            const auto result = to_chars(buf, buf + BufSize, elem, args...);

            global_dummy += static_cast<unsigned int>(result.ptr - buf);
            global_dummy += static_cast<unsigned int>(buf[0]);
        }
    }
    const auto finish = steady_clock::now();

    printf("%6.1f ns | %s\n", duration<double, nano>{finish - start}.count() / (N * K), str);

    for (const auto& elem : vec) {
        const auto result = to_chars(buf, buf + BufSize, elem, args...);
        if (result.ec != errc())
        {
            std::cerr << "To chars failure with: " << elem 
                      << "\n       to_chars value: " << buf << std::endl;
        }

        if constexpr (RT == RoundTrip::Lossy) {
            // skip lossy conversions
        } else {
            T round_trip;
            boost::charconv::from_chars_result from_result;
            if constexpr (std::is_same_v<T, uint32_t> || std::is_same_v<T, uint64_t>)
                from_result = boost::charconv::from_chars(buf, result.ptr, round_trip, args...);
            else
                from_result = boost::charconv::from_chars(buf, result.ptr, round_trip, boost_chars_format_from_RoundTrip(RT));

            if (from_result.ec != errc() || from_result.ptr != result.ptr || round_trip != elem)
            {
                std::cerr << std::setprecision(std::numeric_limits<T>::digits10 + 1)
                          << "Roundtrip failure with: " << elem
                          << "\n          to_chars val: " << buf
                          << "\n        from_chars val: " << round_trip
                          << "\n        from_chars ptr: " << static_cast<std::ptrdiff_t>(from_result.ptr - buf)
                          << "\n          to_chars ptr: " << static_cast<std::ptrdiff_t>(result.ptr - buf) << std::endl;
            }
        }
    }
}

template <typename T>
void test_google_double_conversion_to_string(const char* const str, const vector<T>& values)
{
    using namespace double_conversion;

    const int kBufferSize = 2000;
    char buffer[kBufferSize];
    StringBuilder builder(buffer, kBufferSize);
    const int flags = DoubleToStringConverter::EMIT_POSITIVE_EXPONENT_SIGN;
    DoubleToStringConverter dc(flags, "inf", "nan", 'e', -6, 21, 0, 0, 2);

    vector<string> converted_values(N);

    const auto start = steady_clock::now();
    for (size_t k = 0; k < K; ++k) {
        for (size_t n = 0; n < N; ++n) {
            if constexpr (std::is_same_v<T, double>) {
                dc.ToShortest(values[n], &builder);
            } else {
                dc.ToShortestSingle(values[n], &builder);
            }

            converted_values[n] = builder.Finalize();
            builder.Reset();
        }
    }

    const auto finish = steady_clock::now();

    printf("%6.1f ns | %s\n", duration<double, nano>{finish - start}.count() / (N * K), str);

    for (size_t n = 0; n < N; ++n) {
        T round_trip {};
        auto from_result = from_chars(converted_values[n].data(), converted_values[n].data() + converted_values[n].size(), round_trip);
        verify_fp(round_trip, values[n], 1);
        verify(from_result.ec == errc());
    }
}

template <typename T>
void test_fmt_lib_conversion_to_string(const char* const str, const vector<T>& vec)
{
    char buf[BufSize];

    const auto start = steady_clock::now();
    for (size_t k = 0; k < K; ++k) {
        for (const auto& elem : vec) {
            fmt::format_to(buf, "{}", elem);
            global_dummy += static_cast<unsigned int>(buf[0]);
        }
    }
    const auto finish = steady_clock::now();

    printf("%6.1f ns | %s\n", duration<double, nano>{finish - start}.count() / (N * K), str);

    for (const auto& elem : vec) {
        std::memset(buf, '\0', BufSize);
        fmt::format_to(buf, "{}", elem);

        if constexpr (is_same_v<T, float>) {
            verify_fp(strtof(buf, nullptr), elem);
        } else if constexpr (is_same_v<T, double>) {
            verify_fp(strtod(buf, nullptr), elem);
        } else if constexpr (is_same_v<T, uint32_t>) {
            verify(strtoul(buf, nullptr, 10) == (unsigned long)elem);
        } else if constexpr (is_same_v<T, uint64_t>) {
            verify(strtoull(buf, nullptr, 10) == (unsigned long long)elem);
        }
    }
}

#endif // AVOID_CHARCONV

template <RoundTrip RT, typename T>
vector<char> prepare_strings(const vector<T>& vec) {
    vector<char> output;

    char buf[BufSize];

    for (const auto& elem : vec) {
        int ret;

        if constexpr (RT == RoundTrip::Sci) {
            if constexpr (is_same_v<T, float>) {
                ret = sprintf_wrapper(buf, "%.8e", elem);
            } else {
                ret = sprintf_wrapper(buf, "%.16e", elem);
            }
        } else if constexpr (RT == RoundTrip::Hex) {
            if constexpr (is_same_v<T, float>) {
                ret = sprintf_wrapper(buf, "%.6a", elem);
            } else {
                ret = sprintf_wrapper(buf, "%.13a", elem);
            }
        } else if constexpr (RT == RoundTrip::u32) {
            ret = sprintf_wrapper(buf, "%lu", elem);
        } else {
            static_assert(RT == RoundTrip::u64);
            ret = sprintf_wrapper(buf, "%llu", elem);
        }

        verify(ret != -1);

        output.insert(output.end(), buf, buf + ret + 1); // include null terminator
    }

    return output;
}

template <typename Floating>
void test_strtox(const char* const str, const vector<Floating>& original, const vector<char>& strings) {

    vector<Floating> round_trip(N);

    const auto start = steady_clock::now();
    for (size_t k = 0; k < K; ++k) {
        const char* ptr = strings.data();
        char* endptr    = nullptr;
        for (size_t n = 0; n < N; ++n) {
            if constexpr (is_same_v<Floating, float>) {
                round_trip[n] = strtof(ptr, &endptr);
            } else {
                round_trip[n] = strtod(ptr, &endptr);
            }

            ptr = endptr + 1; // advance past null terminator
        }
    }
    const auto finish = steady_clock::now();

    printf("%6.1f ns | %s\n", duration<double, nano>{finish - start}.count() / (N * K), str);

    verify(round_trip == original);
}

template <typename Floating>
void test_google_double_conversion_from_chars(const char* const str, const vector<Floating>& original, const vector<char>& strings) {

    using namespace double_conversion;

    const char* const last = strings.data() + strings.size();

    vector<Floating> round_trip(N);

    const int flags = StringToDoubleConverter::ALLOW_CASE_INSENSITIVITY;
    int processed;

    auto converter = StringToDoubleConverter(flags, 0.0, 0.0, "inf", "nan");

    const auto start = steady_clock::now();
    for (size_t k = 0; k < K; ++k) {
        const char* ptr = strings.data();
        size_t i = 0;
        while (ptr != last)
        {
            const auto len = strlen(ptr);
            round_trip[i] = converter.StringToDouble(ptr, len, &processed);
            ptr += processed + 1;
            ++i;
        }
    }
    const auto finish = steady_clock::now();

    printf("%6.1f ns | %s\n", duration<double, nano>{finish - start}.count() / (N * K), str);

    for (size_t n = 0; n < N; ++n)
    {
        verify_fp(original[n], round_trip[n], 1);
    }
}

#ifndef AVOID_CHARCONV
vector<char> erase_0x(const vector<char>& strings) {
    vector<char> output;
    output.reserve(strings.size() - 2 * N);

    for (auto i = strings.begin(); i != strings.end();) {
        if (*i == '-') {
            output.push_back('-');
            i += 3; // advance past "-0x";
        } else {
            i += 2; // advance past "0x";
        }

        for (;;) {
            const char c = *i++;
            output.push_back(c);
            if (c == '\0') {
                break;
            }
        }
    }

    return output;
}

template <RoundTrip RT, typename Floating>
void test_from_chars(const char* const str, const vector<Floating>& original, const vector<char>& strings) {

    const char* const last = strings.data() + strings.size();

    vector<Floating> round_trip(N);

    const auto start = steady_clock::now();
    for (size_t k = 0; k < K; ++k) {
        const char* first = strings.data();
        for (size_t n = 0; n < N; ++n) {
            const auto from_result = from_chars(first, last, round_trip[n], chars_format_from_RoundTrip(RT));
            first                  = from_result.ptr + 1; // advance past null terminator
        }
    }
    const auto finish = steady_clock::now();

    printf("%6.1f ns | %s\n", duration<double, nano>{finish - start}.count() / (N * K), str);

    verify(round_trip == original);
}

template <RoundTrip RT, typename Floating>
void test_boost_from_chars(const char* const str, const vector<Floating>& original, const vector<char>& strings) {

    const char* const last = strings.data() + strings.size();

    vector<Floating> round_trip(N);

    const auto start = steady_clock::now();
    for (size_t k = 0; k < K; ++k) {
        const char* first = strings.data();
        for (size_t n = 0; n < N; ++n) {
            const auto from_result = boost::charconv::from_chars(first, last, round_trip[n], boost_chars_format_from_RoundTrip(RT));
            first                  = from_result.ptr + 1; // advance past null terminator
        }
    }
    const auto finish = steady_clock::now();

    printf("%6.1f ns | %s\n", duration<double, nano>{finish - start}.count() / (N * K), str);

    verify(round_trip == original);
}

template <int base, typename Integer>
void test_strtox_integer(const char* const str, const vector<Integer>& original, const vector<char>& strings) {

    vector<Integer> round_trip(N);

    const auto start = steady_clock::now();

    for (size_t k = 0; k < K; ++k) {
        const char* first = strings.data();
        char* end = nullptr;
        for (size_t n = 0; n < N; ++n) {
            if constexpr (std::is_same_v<Integer, uint32_t>) {
                round_trip[n] = std::strtoul(first, &end, base);
            }
            else if constexpr (std::is_same_v<Integer, uint64_t>) {
                round_trip[n] = std::strtoull(first, &end, base);
            }

            first = end + 1; // Advance past the null terminator
        }
    }

    const auto finish = steady_clock::now();

    printf("%6.1f ns | %s\n", duration<double, nano>{finish - start}.count() / (N * K), str);

    verify(round_trip == original);
}

template <int base, typename Integer>
void test_from_chars_integer(const char* const str, [[maybe_unused]] const vector<Integer>& original, const vector<char>& strings) {

    const char* const last = strings.data() + strings.size();

    vector<Integer> round_trip(N);

    const auto start = steady_clock::now();

    for (size_t k = 0; k < K; ++k) {
        const char* first = strings.data();
        for (size_t n = 0; n < N; ++n) {
            const auto from_result = from_chars(first, last, round_trip[n], base);
            first = from_result.ptr + 1; // Advance past the null terminator
        }
    }

    const auto finish = steady_clock::now();

    printf("%6.1f ns | %s\n", duration<double, nano>{finish - start}.count() / (N * K), str);

    verify(round_trip == original);
}

template <int base, typename Integer>
void test_boost_from_chars_integer(const char* const str, [[maybe_unused]] const vector<Integer>& original, const vector<char>& strings) {
    
    const char* const last = strings.data() + strings.size();

    vector<Integer> round_trip(N);

    const auto start = steady_clock::now();
    for (size_t k = 0; k < K; ++k) {
        const char* first = strings.data();
        for (size_t n = 0; n < N; ++n) {
            const auto from_result = boost::charconv::from_chars(first, last, round_trip[n], base);
            first                  = from_result.ptr + 1; // advance past null terminator
        }
    }
    const auto finish = steady_clock::now();

    printf("%6.1f ns | %s\n", duration<double, nano>{finish - start}.count() / (N * K), str);

    verify(round_trip == original);
}

template <typename T, typename Iterator>
bool parse_numbers(Iterator first, Iterator last, std::vector<T>& v)
{
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    namespace phoenix = boost::phoenix;
    using qi::double_;
    using qi::float_;
    using qi::phrase_parse;
    using qi::_1;
    using ascii::space;
    using phoenix::push_back;
    using qi::ulong_;
    using qi::ulong_long;
    using boost::spirit::qi::uint_parser;
    using boost::spirit::qi::real_parser;
    using boost::spirit::qi::real_policies;

    bool r = false;
    
    if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>)
    {
        real_parser<T, qi::strict_real_policies<T>> float_parser;

        if (qi::parse(first, last, float_parser))
        {
            T n;
            auto iter = first;
            size_t i = 0;

            while (qi::parse(iter, last, '\0') && qi::parse(iter, last, float_parser, n))
            {
                v[i] = n;
                ++i;
                first = iter + 1; // Skip null terminator
            }

            return true;
        }
        return false;
    }
    else
    {
        uint_parser<uint64_t, 10, 1, 20> max_uint;

        if (qi::parse(first, last, max_uint))
        {
            uint64_t n;
            auto iter = first;
            size_t i = 0;

            while (qi::parse(iter, last, '\0') && qi::parse(iter, last, max_uint, n))
            {
                v[i] = n;
                ++i;
                first = iter + 1; // Skip null terminator
            }

            return true;
        }
        return false;
    }

    return r;
}

template <typename T>
void test_boost_spirit_qi(const char* const str, BOOST_ATTRIBUTE_UNUSED const vector<T>& original, const vector<char>& strings) {

    const char* const last = strings.data() + strings.size();

    vector<T> round_trip(N);

    const auto start = steady_clock::now();
    for (size_t k = 0; k < K; ++k) {
        const char* first = strings.data();
        parse_numbers(first, last, round_trip);
    }
    const auto finish = steady_clock::now();

    printf("%6.1f ns | %s\n", duration<double, nano>{finish - start}.count() / (N * K), str);

    if constexpr (std::is_same_v<T, uint64_t> || std::is_same_v<T, uint32_t>) {
        for (size_t n = 1; n < N; ++n) {
            verify(original[n] == round_trip[n - 1]);
        }
    } else {
        for (size_t n = 1; n < N; ++n) {
            verify_fp(original[n], round_trip[n - 1], 2);
        }
    }
}

template <typename T>
void test_boost_lexical_cast_parse(const char* const str, const vector<T>& original, const vector<char>& strings) {

    //const char* const last = strings.data() + strings.size();

    vector<T> round_trip(N);

    const auto start = steady_clock::now();
    for (size_t k = 0; k < K; ++k)
    {
        const char* first = strings.data();
        for (size_t n = 0; n < N; ++n)
        {
            const auto len = strlen(first);
            if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>)
                round_trip[n] = boost::lexical_cast<T>(first, len);
            else if constexpr (std::is_same_v<T, uint32_t>)
                round_trip[n] = boost::lexical_cast<unsigned long>(first, len);
            else if constexpr (std::is_same_v<T, uint64_t>)
                round_trip[n] = boost::lexical_cast<unsigned long long>(first, len);
            first += len + 1;
        }
    }
    const auto finish = steady_clock::now();

    printf("%6.1f ns | %s\n", duration<double, nano>{finish - start}.count() / (N * K), str);

    if constexpr (std::is_same_v<T, uint64_t>) {
        verify(round_trip == original);
    } else if constexpr (std::is_floating_point_v<T>) {
        for (size_t n = 0; n < N; ++n) {
            verify_fp(round_trip[n], original[n], 2);
        }
    }
}

/*
template <RoundTrip RT, typename Floating>
void test_boost_from_chars_parser(const char* const str, const vector<Floating>&, const vector<char>& strings) {

    const char* const last = strings.data() + strings.size();

    vector<Floating> round_trip(N);

    bool sign {};
    std::uint64_t significand {};
    std::int64_t  exponent {};

    const auto start = steady_clock::now();
    for (size_t k = 0; k < K; ++k) {
        const char* first = strings.data();
        for (size_t n = 0; n < N; ++n) {
            const auto from_result = boost::charconv::detail::parser(first, last, sign, significand, exponent, boost_chars_format_from_RoundTrip(RT));
            first                  = from_result.ptr + 1; // advance past null terminator
        }
    }
    const auto finish = steady_clock::now();

    printf("%6.1f ns | %s\n", duration<double, nano>{finish - start}.count() / (N * K), str);

    // verify(round_trip == original);
}
*/

#endif // AVOID_CHARCONV

void test_all() {
#if defined(__clang__) && defined(_M_IX86)
    const char* const toolset = "Clang/LLVM x86 + MSVC STL";
#elif defined(__clang__) && defined(_M_X64)
    const char* const toolset = "Clang/LLVM x64 + MSVC STL";
#elif !defined(__clang__) && defined(_M_IX86)
    const char* const toolset = "C1XX/C2 x86 + MSVC STL";
#elif !defined(__clang__) && defined(_M_X64)
    const char* const toolset = "C1XX/C2 x64 + MSVC STL";
#elif defined(__clang__) && defined(__APPLE__) && defined(__arm64__)
    const char* const toolset = "Clang/LLVM Apple ARM + libstdc++";
#elif !defined(__clang__) && defined(__APPLE__) && defined(__arm64__)
    const char* const toolset = "GCC " BOOST_STRINGIZE(__GNUC__) " Apple ARM + libstdc++";
#elif defined(__clang__) && defined(__x86_64__)
    const char* const toolset = "Clang/LLVM " BOOST_STRINGIZE(__clang_major__) " x64 + libstdc++"
#elif defined(__GNUC__) && defined(__x86_64__)
    const char* const toolset = "GCC " BOOST_STRINGIZE(__GNUC__) " x64 + libstdc++" ;
#else
    const char* const toolset = "Unknown Toolset";
#endif
    puts(toolset);


    vector<float> vec_flt;
    vector<double> vec_dbl;
    vector<uint32_t> vec_u32;
    vector<uint64_t> vec_u64;

    {
        mt19937_64 mt64;

        vec_flt.reserve(N);
        while (vec_flt.size() < N) {
            const uint32_t val         = static_cast<uint32_t>(mt64());
            constexpr uint32_t inf_nan = 0x7F800000U;
            if ((val & inf_nan) == inf_nan) {
                continue; // skip INF/NAN
            }
            float flt;
            static_assert(sizeof(flt) == sizeof(val));
            memcpy(&flt, &val, sizeof(flt));
            vec_flt.push_back(flt);
        }

        vec_dbl.reserve(N);
        while (vec_dbl.size() < N) {
            const uint64_t val         = mt64();
            constexpr uint64_t inf_nan = 0x7FF0000000000000ULL;
            if ((val & inf_nan) == inf_nan) {
                continue; // skip INF/NAN
            }
            double dbl;
            static_assert(sizeof(dbl) == sizeof(val));
            memcpy(&dbl, &val, sizeof(dbl));
            vec_dbl.push_back(dbl);
        }

        vec_u32.reserve(N);
        while (vec_u32.size() < N) {
            vec_u32.emplace_back(static_cast<uint32_t>(mt64()));
        }

        vec_u64.reserve(N);
        while (vec_u64.size() < N) {
            vec_u64.emplace_back(mt64());
        }
    }

    puts("\n-----to_chars float-----");

    test_sprintf<RoundTrip::Lossy>("std::sprintf float plain shortest", vec_flt, "%.g");
    test_sprintf<RoundTrip::Lossy>("std::sprintf double plain shortest", vec_dbl, "%.g");

    test_lexical_cast("Boost.lexical_cast float", vec_flt);
    test_lexical_cast("Boost.lexical_cast double", vec_dbl);

    test_boost_spirit_karma("Boost.spirit.karma float", vec_flt);
    test_boost_spirit_karma("Boost.spirit.karma double", vec_dbl);

    /*
    test_sprintf<RoundTrip::Sci>("std::sprintf float scientific 8", vec_flt, "%.8e");
    test_sprintf<RoundTrip::Sci>("std::sprintf double scientific 16", vec_dbl, "%.16e");

    test_sprintf<RoundTrip::Lossy>("std::sprintf float fixed 6 (lossy)", vec_flt, "%f");
    test_sprintf<RoundTrip::Lossy>("std::sprintf double fixed 6 (lossy)", vec_dbl, "%f");

    test_sprintf<RoundTrip::Gen>("std::sprintf float general 9", vec_flt, "%.9g");
    test_sprintf<RoundTrip::Gen>("std::sprintf double general 17", vec_dbl, "%.17g");

    test_sprintf<RoundTrip::Hex>("std::sprintf float hex 6", vec_flt, "%.6a");
    test_sprintf<RoundTrip::Hex>("std::sprintf double hex 13", vec_dbl, "%.13a");
    */
    #ifndef AVOID_CHARCONV
    test_STL_to_chars<RoundTrip::Gen>("std::to_chars float plain shortest", vec_flt);
    test_STL_to_chars<RoundTrip::Gen>("std::to_chars double plain shortest", vec_dbl);
    test_boost_to_chars<RoundTrip::Gen>("Boost.Charconv::to_chars float plain shortest", vec_flt);
    test_boost_to_chars<RoundTrip::Gen>("Boost.Charconv::to_chars double plain shortest", vec_dbl);

    test_google_double_conversion_to_string("double-conversion float plain shortest", vec_flt);
    test_google_double_conversion_to_string("double-conversion double plain shortest", vec_dbl);

    test_fmt_lib_conversion_to_string("{fmt} float scientific", vec_flt);
    test_fmt_lib_conversion_to_string("{fmt} double scientific", vec_dbl);
/*
    test_STL_to_chars<RoundTrip::Sci>("std::to_chars float scientific shortest", vec_flt, chars_format::scientific);
    test_STL_to_chars<RoundTrip::Sci>("std::to_chars double scientific shortest", vec_dbl, chars_format::scientific);
    test_boost_to_chars<RoundTrip::Gen>("Boost.Charconv::to_chars float scientific shortest", vec_flt, boost::charconv::chars_format::scientific);
    test_boost_to_chars<RoundTrip::Gen>("Boost.Charconv::to_chars double scientific shortest", vec_dbl, boost::charconv::chars_format::scientific);

    test_STL_to_chars<RoundTrip::Fix>("std::to_chars float fixed shortest", vec_flt, chars_format::fixed);
    test_STL_to_chars<RoundTrip::Fix>("std::to_chars double fixed shortest", vec_dbl, chars_format::fixed);
    test_boost_to_chars<RoundTrip::Gen>("Boost.Charconv::to_chars float fixed shortest", vec_flt, boost::charconv::chars_format::fixed);
    test_boost_to_chars<RoundTrip::Gen>("Boost.Charconv::to_chars double fixed shortest", vec_dbl, boost::charconv::chars_format::fixed);

    test_STL_to_chars<RoundTrip::Gen>("std::to_chars float general shortest", vec_flt, chars_format::general);
    test_STL_to_chars<RoundTrip::Gen>("std::to_chars double general shortest", vec_dbl, chars_format::general);
    test_boost_to_chars<RoundTrip::Gen>("Boost.Charconv::to_chars float general shortest", vec_flt, boost::charconv::chars_format::general);
    test_boost_to_chars<RoundTrip::Gen>("Boost.Charconv::to_chars double general shortest", vec_dbl, boost::charconv::chars_format::general);

    test_STL_to_chars<RoundTrip::Hex>("std::to_chars float hex shortest", vec_flt, chars_format::hex);
    test_STL_to_chars<RoundTrip::Hex>("std::to_chars double hex shortest", vec_dbl, chars_format::hex);
    test_boost_to_chars<RoundTrip::Hex>("Boost.Charconv::to_chars float hex shortest", vec_flt, boost::charconv::chars_format::hex);
    test_boost_to_chars<RoundTrip::Hex>("Boost.Charconv::to_chars double hex shortest", vec_dbl, boost::charconv::chars_format::hex);

    test_STL_to_chars<RoundTrip::Sci>("std::to_chars float scientific 8", vec_flt, chars_format::scientific, 8);
    test_STL_to_chars<RoundTrip::Sci>("std::to_chars double scientific 16", vec_dbl, chars_format::scientific, 16);
    test_boost_to_chars<RoundTrip::Sci>("Boost.Charconv::to_chars float scientific 8", vec_flt, boost::charconv::chars_format::scientific, 8);
    test_boost_to_chars<RoundTrip::Sci>("Boost.Charconv::to_chars double scientific 16", vec_dbl, boost::charconv::chars_format::scientific, 16);

    test_STL_to_chars<RoundTrip::Lossy>("std::to_chars float fixed 6 (lossy)", vec_flt, chars_format::fixed, 6);
    test_STL_to_chars<RoundTrip::Lossy>("std::to_chars double fixed 6 (lossy)", vec_dbl, chars_format::fixed, 6);
    test_boost_to_chars<RoundTrip::Lossy>("Boost.Charconv::to_chars float fixed 6 (lossy)", vec_flt, boost::charconv::chars_format::fixed, 6);
    test_boost_to_chars<RoundTrip::Lossy>("Boost.Charconv::to_chars double fixed 6 (lossy)", vec_dbl, boost::charconv::chars_format::fixed, 6);

    test_STL_to_chars<RoundTrip::Gen>("std::to_chars float general 9", vec_flt, chars_format::general, 9);
    test_STL_to_chars<RoundTrip::Gen>("std::to_chars double general 17", vec_dbl, chars_format::general, 17);
    test_boost_to_chars<RoundTrip::Gen>("Boost.Charconv::to_chars float general 9", vec_flt, chars_format::general, 9);
    test_boost_to_chars<RoundTrip::Gen>("Boost.Charconv::to_chars double general 17", vec_dbl, chars_format::general, 17);

    test_STL_to_chars<RoundTrip::Hex>("std::to_chars float hex 6", vec_flt, chars_format::hex, 6);
    test_STL_to_chars<RoundTrip::Hex>("std::to_chars double hex 13", vec_dbl, chars_format::hex, 13);
    test_boost_to_chars<RoundTrip::Hex>("Boost.Charconv::to_chars float hex 6", vec_flt, chars_format::hex, 6);
    test_boost_to_chars<RoundTrip::Hex>("Boost.Charconv::to_chars double hex 13", vec_dbl, chars_format::hex, 13);
*/
    #endif // AVOID_CHARCONV

    puts("\n------to_chars int------");

    test_sprintf<RoundTrip::Gen>("std::sprintf uint32_t", vec_u32, "%lu");
    test_sprintf<RoundTrip::Gen>("std::sprintf uint64_t", vec_u64, "%llu");

    test_lexical_cast("Boost.lexical_cast uint32_t", vec_u32);
    test_lexical_cast("Boost.lexical_cast uint64_t", vec_u64);

    test_boost_spirit_karma("Boost.spirit.karma uint32_t", vec_u32);
    test_boost_spirit_karma("Boost.spirit.karma uint64_t", vec_u64);

    test_STL_to_chars<RoundTrip::Gen>("std::to_chars uint32_t", vec_u32, 10);
    test_STL_to_chars<RoundTrip::Gen>("std::to_chars uint64_t", vec_u64, 10);

    test_boost_to_chars<RoundTrip::Gen>("Boost.Charconv::to_chars uint32_t", vec_u32, 10);
    test_boost_to_chars<RoundTrip::Gen>("Boost.Charconv::to_chars uint64_t", vec_u64, 10);

    test_fmt_lib_conversion_to_string("{fmt} uint32_t", vec_u32);
    test_fmt_lib_conversion_to_string("{fmt} uint64_t", vec_u64);

    puts("\n----from_chars float----");

    const vector<char> strings_sci_flt = prepare_strings<RoundTrip::Sci>(vec_flt);
    const vector<char> strings_sci_dbl = prepare_strings<RoundTrip::Sci>(vec_dbl);

    const vector<char> strings_hex_flt = prepare_strings<RoundTrip::Hex>(vec_flt);
    const vector<char> strings_hex_dbl = prepare_strings<RoundTrip::Hex>(vec_dbl);

    const vector<char> strings_u32 = prepare_strings<RoundTrip::u32>(vec_u32);
    const vector<char> strings_u64 = prepare_strings<RoundTrip::u64>(vec_u64);

    test_strtox("std::strtof float scientific", vec_flt, strings_sci_flt);
    test_strtox("std::strtod double scientific", vec_dbl, strings_sci_dbl);

    test_boost_lexical_cast_parse("Boost.lexical_cast float scientific", vec_flt, strings_sci_flt);
    test_boost_lexical_cast_parse("Boost.lexical_cast double scientific", vec_dbl, strings_sci_dbl);

    test_boost_spirit_qi<float>("Boost.Spirit.Qi float scientific", vec_flt, strings_sci_flt);
    test_boost_spirit_qi<double>("Boost.Spirit.Qi double scientific", vec_dbl, strings_sci_dbl);

    test_from_chars<RoundTrip::Sci>("std::from_chars float scientific", vec_flt, strings_sci_flt);
    test_from_chars<RoundTrip::Sci>("std::from_chars double scientific", vec_dbl, strings_sci_dbl);

    test_boost_from_chars<RoundTrip::Sci>("Boost.Charconv::from_chars float scientific", vec_flt, strings_sci_flt);
    test_boost_from_chars<RoundTrip::Sci>("Boost.Charconv::from_chars double scientific", vec_dbl, strings_sci_dbl);

    test_google_double_conversion_from_chars("double-conversion float scientific", vec_flt, strings_sci_flt);
    test_google_double_conversion_from_chars("double-conversion double scientific", vec_dbl, strings_sci_dbl);

    //test_strtox("std::strtof float hex", vec_flt, strings_hex_flt);
    //test_strtox("std::strtod double hex", vec_dbl, strings_hex_dbl);

    //test_boost_from_chars_parser<RoundTrip::Sci>("Boost.Charconv::from_chars::parser float scientific", vec_flt, strings_sci_flt);
    //test_boost_from_chars_parser<RoundTrip::Sci>("Boost.Charconv::from_chars::parser double scientific", vec_dbl, strings_sci_dbl);

    //test_from_chars<RoundTrip::Hex>("std::from_chars float hex", vec_flt, erase_0x(strings_hex_flt));
    //test_from_chars<RoundTrip::Hex>("std::from_chars double hex", vec_dbl, erase_0x(strings_hex_dbl));
    //test_boost_from_chars<RoundTrip::Hex>("Boost.Charconv::from_chars float hex", vec_flt, erase_0x(strings_hex_flt));
    //test_boost_from_chars<RoundTrip::Hex>("Boost.Charconv::from_chars double hex", vec_dbl, erase_0x(strings_hex_dbl));

    puts("\n-----from_chars int-----");

    test_strtox_integer<10>("std::strtoul uint32_t", vec_u32, strings_u32);
    test_strtox_integer<10>("std::strtoull uint64_t", vec_u64, strings_u64);

    test_boost_lexical_cast_parse("Boost.lexical_cast uint32_t", vec_u32, strings_u32);
    test_boost_lexical_cast_parse("Boost.lexical_cast uint64_t", vec_u64, strings_u64);

    test_boost_spirit_qi("Boost.Spirit.Qi uint32_t", vec_u32, strings_u32);
    test_boost_spirit_qi("Boost.Spirit.Qi uint64_t", vec_u64, strings_u64);

    test_from_chars_integer<10>("std::from_chars uint32_t", vec_u32, strings_u32);
    test_from_chars_integer<10>("std::from_chars uint64_t", vec_u64, strings_u64);

    test_boost_from_chars_integer<10>("Boost.Charconv::from_chars uint32_t", vec_u32, strings_u32);
    test_boost_from_chars_integer<10>("Boost.Charconv::from_chars uint64_t", vec_u64, strings_u64);

    printf("global_dummy: %u\n", global_dummy);
}

int main()
{
    try {
        test_all();
    } catch (const exception& e) {
        printf("Exception: %s\n", e.what());
    } catch (...) {
        printf("Unknown exception.\n");
    }

    return 1;
}

#else

int main()
{
    std::cerr << "Benchmarks not run" << std::endl;
    return 1;
}

#endif
