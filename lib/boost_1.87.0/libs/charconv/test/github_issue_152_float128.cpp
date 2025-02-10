// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/charconv/detail/config.hpp>
#include <ostream>

#ifdef BOOST_CHARCONV_HAS_QUADMATH

std::ostream& operator<<( std::ostream& os, __float128 v )
{
    char buffer[ 256 ] {};
    quadmath_snprintf(buffer, sizeof(buffer), "%Qg", v);
    os << buffer;
    return os;
}

#ifdef BOOST_CHARCONV_HAS_STDFLOAT128
std::ostream& operator<<( std::ostream& os, std::float128_t v )
{
    os << static_cast<__float128>(v);
    return os;
}
#endif

#include <boost/charconv.hpp>
#include <boost/core/lightweight_test.hpp>
#include <system_error>
#include <limits>
#include <random>
#include <array>
#include <cstdint>

constexpr std::size_t N = 1024;
static std::mt19937_64 rng(42);

template <typename T>
void test_non_finite()
{
    constexpr std::array<T, 4> values = {{INFINITY, -INFINITY, NAN, -NAN}};

    for (const auto val : values)
    {
        char buffer[2];
        auto r = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), val);
        BOOST_TEST(r.ec == std::errc::value_too_large);
    }

    char inf_buffer[3];
    auto r_inf = boost::charconv::to_chars(inf_buffer, inf_buffer + 3, values[0]);
    BOOST_TEST(r_inf);
    BOOST_TEST(!std::memcmp(inf_buffer, "inf", 3));

    char nan_buffer[3];
    auto r_nan = boost::charconv::to_chars(nan_buffer, nan_buffer + 3, values[2]);
    BOOST_TEST(r_nan);
    BOOST_TEST(!std::memcmp(nan_buffer, "nan", 3));

    char neg_nan_buffer[9];
    auto r_neg_nan = boost::charconv::to_chars(neg_nan_buffer, neg_nan_buffer + 9, values[3]);
    BOOST_TEST(r_neg_nan);
    BOOST_TEST(!std::memcmp(neg_nan_buffer, "-nan(ind)", 9));
}

template <typename T>
void test_min_buffer_size()
{
    // No guarantees are made for fixed, especially in this domain
    auto formats = {boost::charconv::chars_format::hex,
                    boost::charconv::chars_format::scientific,
                    boost::charconv::chars_format::general};

    for (const auto format : formats)
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            char buffer[boost::charconv::limits<T>::max_chars10];
            const auto value = static_cast<T>(rng());

            auto r = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), value, format);
            if (!BOOST_TEST(r))
            {
                std::cerr << "Overflow for: " << value << std::endl; // LCOV_EXCL_LINE
            }
        }
    }
}

int main()
{
    test_non_finite<__float128>();
    test_min_buffer_size<__float128>();

    #ifdef BOOST_CHARCONV_HAS_STDFLOAT128
    test_non_finite<std::float128_t>();
    test_min_buffer_size<std::float128_t>();
    #endif

    return boost::report_errors();
}

#else

int main()
{
    return 0;
}

#endif
