// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/charconv.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <system_error>
#include <limits>
#include <random>
#include <array>
#include <cstdint>
#include <iomanip>

constexpr std::size_t N = 1024;
static std::mt19937_64 rng(42);

template <typename T>
void test_non_finite()
{
    constexpr std::array<T, 6> values = {{std::numeric_limits<T>::infinity(), -std::numeric_limits<T>::infinity(),
                                         std::numeric_limits<T>::quiet_NaN(), -std::numeric_limits<T>::quiet_NaN(),
                                         std::numeric_limits<T>::signaling_NaN(), -std::numeric_limits<T>::signaling_NaN()}};

    for (const auto val : values)
    {
        char buffer[2];
        auto r = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), val);
        BOOST_TEST(r.ec == std::errc::value_too_large);
    }

    char inf_buffer[3];
    constexpr T inf_val = std::numeric_limits<T>::infinity();
    auto r_inf = boost::charconv::to_chars(inf_buffer, inf_buffer + 3, inf_val);
    BOOST_TEST(r_inf);
    BOOST_TEST(!std::memcmp(inf_buffer, "inf", 3));

    char nan_buffer[3];
    constexpr T nan_val = std::numeric_limits<T>::quiet_NaN();
    auto r_nan = boost::charconv::to_chars(nan_buffer, nan_buffer + 3, nan_val);
    BOOST_TEST(r_nan);
    BOOST_TEST(!std::memcmp(nan_buffer, "nan", 3));

    char neg_nan_buffer[9];
    auto r_neg_nan = boost::charconv::to_chars(neg_nan_buffer, neg_nan_buffer + 9, -nan_val);
    BOOST_TEST(r_neg_nan);
    BOOST_TEST(!std::memcmp(neg_nan_buffer, "-nan(ind)", 9));

    char snan_buffer[9];
    constexpr T snan_val = std::numeric_limits<T>::signaling_NaN();
    auto r_snan = boost::charconv::to_chars(snan_buffer, snan_buffer + 9, snan_val);
    BOOST_TEST(r_snan);
    BOOST_TEST(!std::memcmp(snan_buffer, "nan(snan)", 9));
}

template <typename T>
void test_non_finite_fixed_precision()
{
    constexpr std::array<T, 6> values = {{std::numeric_limits<T>::infinity(), -std::numeric_limits<T>::infinity(),
                                         std::numeric_limits<T>::quiet_NaN(), -std::numeric_limits<T>::quiet_NaN(),
                                         std::numeric_limits<T>::signaling_NaN(), -std::numeric_limits<T>::signaling_NaN()}};

    const auto formats = {boost::charconv::chars_format::scientific, boost::charconv::chars_format::hex,
                              boost::charconv::chars_format::general, boost::charconv::chars_format::fixed};

    for (auto format : formats)
    {
        for (const auto val : values)
        {
            char buffer[2];
            auto r = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), val, format, 5);
            BOOST_TEST(r.ec == std::errc::value_too_large);
        }

        char inf_buffer[3];
        constexpr T inf_val = std::numeric_limits<T>::infinity();
        auto r_inf = boost::charconv::to_chars(inf_buffer, inf_buffer + 3, inf_val, format, 5);
        BOOST_TEST(r_inf);
        BOOST_TEST(!std::memcmp(inf_buffer, "inf", 3));

        char nan_buffer[3];
        constexpr T nan_val = std::numeric_limits<T>::quiet_NaN();
        auto r_nan = boost::charconv::to_chars(nan_buffer, nan_buffer + 3, nan_val, format, 5);
        BOOST_TEST(r_nan);
        BOOST_TEST(!std::memcmp(nan_buffer, "nan", 3));

        char neg_nan_buffer[9];
        auto r_neg_nan = boost::charconv::to_chars(neg_nan_buffer, neg_nan_buffer + 9, -nan_val, format, 5);
        BOOST_TEST(r_neg_nan);
        BOOST_TEST(!std::memcmp(neg_nan_buffer, "-nan(ind)", 9));

        char snan_buffer[9];
        constexpr T snan_val = std::numeric_limits<T>::signaling_NaN();
        auto r_snan = boost::charconv::to_chars(snan_buffer, snan_buffer + 9, snan_val, format, 5);
        BOOST_TEST(r_snan);
        BOOST_TEST(!std::memcmp(snan_buffer, "nan(snan)", 9));
    }
};

#ifdef BOOST_MSVC
# pragma warning(push)
# pragma warning(disable: 4127)
#endif

template <typename T>
void test_min_buffer_size()
{
    #if defined(_WIN32)
    std::uniform_real_distribution<T> dist((std::numeric_limits<T>::min)(), (std::numeric_limits<T>::max)());
    #else
    std::uniform_real_distribution<T> dist((std::numeric_limits<T>::lowest)(), (std::numeric_limits<T>::max)());
    #endif

    // No guarantees are made for fixed, especially in this domain
    auto formats = {boost::charconv::chars_format::hex,
                                         boost::charconv::chars_format::scientific,
                                         boost::charconv::chars_format::general};

    int format_int = 0;
    for (const auto format : formats)
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            char buffer[boost::charconv::limits<T>::max_chars10];
            const T value = dist(rng);

            if (!std::isnormal(value))
            {
                continue;
            }

            auto r = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), value, format);
            if (!BOOST_TEST(r))
            {
                // LCOV_EXCL_START
                std::cerr << std::setprecision(std::numeric_limits<T>::max_digits10) << "Overflow for: " << value
                          << "\nFormat: " << format_int
                          << "\nBuffer size: " << sizeof(buffer) << std::endl;
                // LCOV_EXCL_STOP
            }
        }
        ++format_int;
    }
}

#ifdef BOOST_MSVC
# pragma warning(pop)
#endif

#if BOOST_CHARCONV_LDBL_BITS > 64
void test_failed_values()
{
    // No guarantees are made for fixed, especially in this domain
    // TODO(mborland): Add hex once https://github.com/boostorg/charconv/issues/154 is fixed
    auto formats = {boost::charconv::chars_format::scientific,
                                         boost::charconv::chars_format::general};

    std::array<long double, 2> failed_values = {{6.93880126833169422964e+4931L, 9.14517491001980558957e+4931L}};

    int format_int = 0;
    for (const auto format : formats)
    {
        for (const auto value : failed_values)
        {
            char buffer[boost::charconv::limits<long double>::max_chars10];

            if (!std::isnormal(value))
            {
                continue;
            }

            auto r = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), value, format);
            if (!BOOST_TEST(r))
            {
                // LCOV_EXCL_START
                std::cerr << std::setprecision(std::numeric_limits<long double>::max_digits10) << "Overflow for: " << value
                          << "\nFormat: " << format_int
                          << "\nBuffer size: " << sizeof(buffer) << std::endl;
                // LCOV_EXCL_STOP
            }
        }
        ++format_int;
    }
}
#endif

int main()
{
    test_non_finite<float>();
    test_non_finite<double>();
    #ifdef BOOST_CHARCONV_HAS_FLOAT16
    test_non_finite<std::float16_t>();
    #endif
    #ifdef BOOST_CHARCONV_HAS_FLOAT32
    test_non_finite<std::float32_t>();
    #endif
    #ifdef BOOST_CHARCONV_HAS_FLOAT64
    test_non_finite<std::float64_t>();
    #endif
    #ifdef BOOST_CHARCONV_HAS_BRAINFLOAT16
    test_non_finite<std::bfloat16_t>();
    #endif

    test_non_finite_fixed_precision<float>();
    test_non_finite_fixed_precision<double>();
    #ifdef BOOST_CHARCONV_HAS_FLOAT16
    test_non_finite_fixed_precision<std::float16_t>();
    #endif
    #ifdef BOOST_CHARCONV_HAS_FLOAT32
    test_non_finite_fixed_precision<std::float32_t>();
    #endif
    #ifdef BOOST_CHARCONV_HAS_FLOAT64
    test_non_finite_fixed_precision<std::float64_t>();
    #endif
    #ifdef BOOST_CHARCONV_HAS_BRAINFLOAT16
    test_non_finite_fixed_precision<std::bfloat16_t>();
    #endif

    test_min_buffer_size<float>();
    test_min_buffer_size<double>();
    #ifdef BOOST_CHARCONV_HAS_FLOAT32
    test_min_buffer_size<std::float32_t>();
    #endif
    #ifdef BOOST_CHARCONV_HAS_FLOAT64
    test_min_buffer_size<std::float64_t>();
    #endif

    #if BOOST_CHARCONV_LDBL_BITS > 64
    test_failed_values();
    #endif

    #ifndef BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE
    test_non_finite<long double>();
    test_non_finite_fixed_precision<long double>();
    test_min_buffer_size<long double>();
    #endif

    return boost::report_errors();
}
