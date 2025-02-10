// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/charconv.hpp>
#include <boost/core/lightweight_test.hpp>

template <typename T>
void overflow_spot_value(const std::string& buffer, boost::charconv::chars_format fmt = boost::charconv::chars_format::general)
{
    auto v = static_cast<T>(42.L);
    auto r = boost::charconv::from_chars(buffer.c_str(), buffer.c_str() + std::strlen(buffer.c_str()), v, fmt);

    BOOST_TEST(v == static_cast<T>(42.L)) && BOOST_TEST(r.ec == std::errc::result_out_of_range);
}

template <typename T>
void test()
{
    const auto format_list = {boost::charconv::chars_format::general, boost::charconv::chars_format::scientific, boost::charconv::chars_format::hex};

    for (const auto format : format_list)
    {
        if (format != boost::charconv::chars_format::hex)
        {
            overflow_spot_value<T>("1e99999", format);
            overflow_spot_value<T>("-1e99999", format);
            overflow_spot_value<T>("1e-99999", format);
            overflow_spot_value<T>("-1.0e-99999", format);
        }
        else
        {
            overflow_spot_value<T>("1p99999", format);
            overflow_spot_value<T>("-1p99999", format);
            overflow_spot_value<T>("1p-99999", format);
            overflow_spot_value<T>("-1.0p-99999", format);
        }
    }
}

int main()
{
    test<float>();
    test<double>();

    #ifndef BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE
    test<long double>();
    #endif

    #ifdef BOOST_CHARCONV_HAS_FLOAT128
    test<__float128>();
    #endif

    #ifdef BOOST_CHARCONV_HAS_FLOAT16
    test<std::float16_t>();
    #endif
    #ifdef BOOST_CHARCONV_HAS_FLOAT32
    test<std::float32_t>();
    #endif
    #ifdef BOOST_CHARCONV_HAS_FLOAT64
    test<std::float64_t>();
    #endif
    #if defined(BOOST_CHARCONV_HAS_STDFLOAT128) && defined(BOOST_CHARCONV_HAS_FLOAT128)
    test<std::float128_t>();
    #endif
    #ifdef BOOST_CHARCONV_HAS_BRAINFLOAT16
    test<std::bfloat16_t>();
    #endif

    return boost::report_errors();
}
