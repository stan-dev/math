// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/charconv/issues/186

#include <boost/charconv.hpp>
#include <boost/core/lightweight_test.hpp>
#include <iostream>
#include <string>

template <typename T>
void force_overflow(T value)
{
    for (int precision = 0; precision < 1000; ++precision)
    {
        char buffer [20] {};
        //std::cerr << "Current precision: " << precision << std::endl;
        boost::charconv::to_chars(buffer, buffer + sizeof(buffer), value, boost::charconv::chars_format::fixed, precision);
    }
}

template <typename T>
void force_overflow_all_formats(T value)
{
    const auto formats = {boost::charconv::chars_format::general,
                          boost::charconv::chars_format::fixed,
                          boost::charconv::chars_format::scientific,
                          boost::charconv::chars_format::hex};

    for (const auto format : formats)
    {
        char buffer[20]; // Small enough it should encounter overflows

        boost::charconv::to_chars(buffer, buffer + sizeof(buffer), value, format);

        // Also try with precisions
        for (int precision = -1; precision < 1000; ++precision)
        {
            boost::charconv::to_chars(buffer, buffer + sizeof(buffer), value, format, precision);
        }
    }
}

template <typename T>
void spot_check_result(T value, const char* string)
{
    char buffer[21] {};
    const auto r = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, value, boost::charconv::chars_format::fixed, 0);
    BOOST_TEST(r);
    *r.ptr = '\0';

    BOOST_TEST_CSTR_EQ(buffer, string);
}

int main()
{
    // See fuzzing/fuzz_to_chars_float
    // Validation values calculated using GCC 13.2 on godbolt: https://godbolt.org/z/8ns7EreGW
    force_overflow(444444444442.0);
    force_overflow(444444444442.0F);
    spot_check_result(444444444442.0, "444444444442");

    force_overflow(33333333333333333333.0);
    force_overflow(33333333333333333333.0F);
    spot_check_result(33333333333333333333.0, "33333333333333331968");

    force_overflow(66666666666666666666.0);
    force_overflow(66666666666666666666.0F);
    spot_check_result(66666666666666666666.0, "66666666666666663936");

    force_overflow_all_formats(55555555095580.825555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555);

    return boost::report_errors();
}
