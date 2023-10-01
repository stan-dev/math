//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale/conversion.hpp>
#include <boost/locale/generator.hpp>
#include <boost/locale/info.hpp>
#include <boost/locale/localization_backend.hpp>
#include "boostLocale/test/tools.hpp"
#include "boostLocale/test/unit_test.hpp"
#include "case_convert_test.hpp"
#include <iomanip>
#include <iostream>

template<typename CharType>
void test_char()
{
    using boost::locale::case_convert_test::test_one;

    boost::locale::generator gen;

    std::cout << "- Testing at least C" << std::endl;
    std::locale l = gen("en_US.UTF-8");
    test_one<CharType>(l, "Hello World i", "hello world i", "HELLO WORLD I");

    std::string name;

    name = get_std_name("en_US.UTF-8");
    if(!name.empty()) {
        std::cout << "- Testing " << name << std::endl;
        l = gen(name);
        test_one<CharType>(l, "Façade", "façade", "FAÇADE");
    } else
        std::cout << "- en_US.UTF-8 is not supported, skipping" << std::endl; // LCOV_EXCL_LINE

    name = get_std_name("en_US.ISO8859-1");
    if(!name.empty()) {
        std::cout << "Testing " << name << std::endl;
        l = gen(name);
        test_one<CharType>(l, "Hello World", "hello world", "HELLO WORLD");
        test_one<CharType>(l, "Façade", "façade", "FAÇADE");
    } else
        std::cout << "- en_US.ISO8859-1 is not supported, skipping" << std::endl; // LCOV_EXCL_LINE

    std::string real_name;
    name = get_std_name("tr_TR.UTF-8", &real_name);
    if(!name.empty()) {
        std::cout << "Testing " << name << std::endl;
        if(std::use_facet<std::ctype<wchar_t>>(std::locale(real_name)).toupper(L'i') != L'I') {
            l = gen(name);
            test_one<CharType>(l, "i", "i", "İ");
        } else
            std::cout << "Standard library does not support this locale's case conversion correctly" << std::endl;
    } else
        std::cout << "- tr_TR.UTF-8 is not supported, skipping" << std::endl; // LCOV_EXCL_LINE
}

BOOST_LOCALE_DISABLE_UNREACHABLE_CODE_WARNING
void test_main(int /*argc*/, char** /*argv*/)
{
#ifdef BOOST_LOCALE_NO_STD_BACKEND
    std::cout << "STD Backend is not build... Skipping\n";
    return;
#endif
    boost::locale::localization_backend_manager mgr = boost::locale::localization_backend_manager::global();
    mgr.select("std");
    boost::locale::localization_backend_manager::global(mgr);

    std::cout << "Testing char" << std::endl;
    test_char<char>();
    std::cout << "Testing wchar_t" << std::endl;
    test_char<wchar_t>();
#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
    std::cout << "Testing char16_t" << std::endl;
    test_char<char16_t>();
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
    std::cout << "Testing char32_t" << std::endl;
    test_char<char32_t>();
#endif
}

// boostinspect:noascii
