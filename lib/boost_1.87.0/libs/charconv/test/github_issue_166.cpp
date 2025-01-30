// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/cppalliance/charconv/issues/166

#include <boost/charconv.hpp>
#include <boost/core/lightweight_test.hpp>

template <typename T>
void simple_test()
{
    for (int i = -1; i > -5; --i)
    {
        char buffer[256] {};
        boost::charconv::to_chars(buffer, buffer + sizeof(buffer), static_cast<T>(9.99999L), boost::charconv::chars_format::fixed, i);
        BOOST_TEST_CSTR_EQ(buffer, "9.999990");
    }
}

int main()
{
    simple_test<float>();
    simple_test<double>();
    #if BOOST_CHARCONV_LDBL_BITS == 64
    simple_test<long double>();
    #endif

    #ifdef BOOST_CHARCONV_HAS_FLOAT32
    simple_test<std::float32_t>();
    #endif
    #ifdef BOOST_CHARCONV_HAS_FLOAT64
    simple_test<std::float64_t>();
    #endif

    return boost::report_errors();
}
