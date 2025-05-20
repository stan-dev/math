//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale/conversion.hpp>
#include <boost/locale/generator.hpp>
#include <boost/locale/info.hpp>
#include "boostLocale/test/tools.hpp"
#include "boostLocale/test/unit_test.hpp"
#include <iomanip>
#include <iostream>

template<typename Char>
void test_normc(std::basic_string<Char> orig, std::basic_string<Char> normal, boost::locale::norm_type type)
{
    std::locale l = boost::locale::generator().generate("en_US.UTF-8");
    TEST_EQ(normalize(orig, type, l), normal);
    TEST_EQ(normalize(orig.c_str(), type, l), normal);
    TEST_EQ(normalize(orig.c_str(), orig.c_str() + orig.size(), type, l), normal);
}

void test_norm(std::string orig, std::string normal, boost::locale::norm_type type)
{
    test_normc<char>(orig, normal, type);
    test_normc<wchar_t>(to<wchar_t>(orig), to<wchar_t>(normal), type);
#ifndef BOOST_LOCALE_NO_CXX20_STRING8
    test_normc<char8_t>(to<char8_t>(orig), to<char8_t>(normal), type);
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
    test_normc<char16_t>(to<char16_t>(orig), to<char16_t>(normal), type);
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
    test_normc<char32_t>(to<char32_t>(orig), to<char32_t>(normal), type);
#endif
}

#define TEST_A(Chr, how, source, dest)                                                             \
    do {                                                                                           \
        const boost::locale::info& inf = std::use_facet<boost::locale::info>(std::locale());       \
        std::cout << "Testing " #how " for " #Chr ", lang=" << inf.language();                     \
        if(std::string("char") == #Chr)                                                            \
            std::cout << " charset=" << inf.encoding();                                            \
        std::cout << std::endl;                                                                    \
        std::basic_string<Chr> source_s = (source), dest_s = (dest);                               \
        TEST_EQ(boost::locale::how(source_s), dest_s);                                             \
        TEST_EQ(boost::locale::how(source_s.c_str()), dest_s);                                     \
        TEST_EQ(boost::locale::how(source_s.c_str(), source_s.c_str() + source_s.size()), dest_s); \
        BOOST_LOCALE_START_CONST_CONDITION                                                         \
    } while(0) BOOST_LOCALE_END_CONST_CONDITION

#define TEST_ALL_CASES                                    \
    do {                                                  \
        eight_bit = true;                                 \
        std::locale::global(gen("en_US.UTF-8"));          \
        TEST_V(to_upper, "grüßen i", "GRÜSSEN I");        \
        TEST_V(to_lower, "Façade", "façade");             \
        TEST_V(to_title, "façadE world", "Façade World"); \
        TEST_V(fold_case, "Hello World", "hello world");  \
        std::locale::global(gen("tr_TR.UTF-8"));          \
        eight_bit = false;                                \
        TEST_V(to_upper, "i", "İ");                       \
        TEST_V(to_lower, "İ", "i");                       \
        BOOST_LOCALE_START_CONST_CONDITION                \
    } while(0) BOOST_LOCALE_END_CONST_CONDITION

BOOST_LOCALE_DISABLE_UNREACHABLE_CODE_WARNING
void test_main(int /*argc*/, char** /*argv*/)
{
#ifndef BOOST_LOCALE_WITH_ICU
    std::cout << "ICU is not build... Skipping\n";
    return;
#endif
    {
        using namespace boost::locale;
        std::cout << "Testing Unicode normalization" << std::endl;
        test_norm("\xEF\xAC\x81", "\xEF\xAC\x81", norm_nfd); // ligature fi
        test_norm("\xEF\xAC\x81", "\xEF\xAC\x81", norm_nfc);
        test_norm("\xEF\xAC\x81", "fi", norm_nfkd);
        test_norm("\xEF\xAC\x81", "fi", norm_nfkc);
        test_norm("ä", "ä", norm_nfd); // ä to a and accent
        test_norm("ä", "ä", norm_nfc);
    }

    boost::locale::generator gen;
    bool eight_bit = true;

#define TEST_V(how, source_s, dest_s)                                \
    do {                                                             \
        TEST_A(char, how, source_s, dest_s);                         \
        if(eight_bit) {                                              \
            std::locale tmp = std::locale();                         \
            std::locale::global(gen("en_US.ISO8859-1"));             \
            TEST_A(char, how, to<char>(source_s), to<char>(dest_s)); \
            std::locale::global(tmp);                                \
        }                                                            \
        BOOST_LOCALE_START_CONST_CONDITION                           \
    } while(0) BOOST_LOCALE_END_CONST_CONDITION

    TEST_ALL_CASES;
#undef TEST_V

#define TEST_V(how, source_s, dest_s) TEST_A(wchar_t, how, to<wchar_t>(source_s), to<wchar_t>(dest_s))
    TEST_ALL_CASES;
#undef TEST_V

#ifndef BOOST_LOCALE_NO_CXX20_STRING8
#    define TEST_V(how, source_s, dest_s) TEST_A(char8_t, how, to<char8_t>(source_s), to<char8_t>(dest_s))
    TEST_ALL_CASES;
#    undef TEST_V
#endif

#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
#    define TEST_V(how, source_s, dest_s) TEST_A(char16_t, how, to<char16_t>(source_s), to<char16_t>(dest_s))
    TEST_ALL_CASES;
#    undef TEST_V
#endif

#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
#    define TEST_V(how, source_s, dest_s) TEST_A(char32_t, how, to<char32_t>(source_s), to<char32_t>(dest_s))
    TEST_ALL_CASES;
#    undef TEST_V
#endif
}

// boostinspect:noascii
