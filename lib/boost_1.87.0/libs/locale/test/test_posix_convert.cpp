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
#include <wctype.h>

template<typename CharType>
void test_char()
{
    using boost::locale::case_convert_test::test_one;
    boost::locale::generator gen;

    std::cout << "- Testing at least C" << std::endl;
    std::locale l = gen("C.UTF-8");
    test_one<CharType>(l, "Hello World i", "hello world i", "HELLO WORLD I");
    boost::locale::case_convert_test::test_no_op_title_case<CharType>(l, "Hello world i");

    std::string name = "en_US.UTF-8";
    if(has_posix_locale(name)) {
        std::cout << "- Testing " << name << std::endl;
        l = gen(name);
        test_one<CharType>(l, "Façade", "façade", "FAÇADE");
        boost::locale::case_convert_test::test_no_op_title_case<CharType>(l, "Hello world i");
    }

    name = "en_US.ISO8859-1";
    if(has_posix_locale(name)) {
        std::cout << "Testing " << name << std::endl;
        l = gen(name);
        test_one<CharType>(l, "Hello World", "hello world", "HELLO WORLD");
        boost::locale::case_convert_test::test_no_op_title_case<CharType>(l, "Hello world i");
#if defined(__APPLE__) || defined(__FreeBSD__)
        if(sizeof(CharType) != 1)
#endif
            test_one<CharType>(l, "Façade", "façade", "FAÇADE");
    }

    name = "tr_TR.UTF-8";
    if(!has_posix_locale(name))
        std::cout << "- " << name << " is not supported, skipping" << std::endl; // LCOV_EXCL_LINE
    else {
        std::cout << "Testing " << name << std::endl;
        locale_holder cl(newlocale(LC_ALL_MASK, name.c_str(), nullptr));
        TEST_REQUIRE(cl);
#ifndef BOOST_LOCALE_NO_POSIX_BACKEND
        if(towupper_l(L'i', cl) == 0x130)
            test_one<CharType>(gen(name), "i", "i", "İ");
        else
            std::cout << "  Turkish locale is not supported well" << std::endl; // LCOV_EXCL_LINE
#endif
    }
}

BOOST_LOCALE_DISABLE_UNREACHABLE_CODE_WARNING
void test_main(int /*argc*/, char** /*argv*/)
{
#ifdef BOOST_LOCALE_NO_POSIX_BACKEND
    std::cout << "POSIX Backend is not build... Skipping\n";
    return;
#endif
    boost::locale::localization_backend_manager mgr = boost::locale::localization_backend_manager::global();
    mgr.select("posix");
    boost::locale::localization_backend_manager::global(mgr);

    std::cout << "Testing char" << std::endl;
    test_char<char>();
    std::cout << "Testing wchar_t" << std::endl;
    test_char<wchar_t>();
}

// boostinspect:noascii
