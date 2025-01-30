// Copyright 2024 Alexander Grund
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <locale>
#include <iostream>
#include <boost/charconv.hpp>
#include <boost/charconv/detail/fallback_routines.hpp>
#include <boost/core/lightweight_test.hpp>

template <typename T>
void test()
{
    const char buffer[] = "1.1897e+2";
    constexpr auto valdiation_value = static_cast<T>(1.1897e+2L);

    try
    {
        #ifdef BOOST_MSVC
        std::locale::global(std::locale("German"));
        #else
        std::locale::global(std::locale("de_DE.UTF-8"));
        #endif
    }
    // LCOV_EXCL_START
    catch (...)
    {
        std::cerr << "Locale not installed. Skipping test." << std::endl;
        return;
    }
    // LCOV_EXCL_STOP

    T v = 0;
    auto r = boost::charconv::detail::from_chars_strtod(buffer, buffer + sizeof(buffer), v);
    BOOST_TEST(r);

    std::locale::global(std::locale::classic());
    T v2 = 0;
    auto r2 = boost::charconv::detail::from_chars_strtod(buffer, buffer + sizeof(buffer), v2);
    BOOST_TEST(r2);

    BOOST_TEST_EQ(v, v2);
    BOOST_TEST_EQ(v, valdiation_value);
    BOOST_TEST(r.ptr == r2.ptr);
}

// See: https://stackoverflow.com/questions/1745045/stdlocale-breakage-on-macos-10-6-with-lang-en-us-utf-8
#if !defined(__APPLE__) || (defined(__APPLE__) && defined(__clang__))

int main()
{
    test<float>();
    test<double>();

    #ifndef BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE
    test<long double>();
    #endif

    return boost::report_errors();
}

#else

int main()
{
    return 0;
}

#endif
