//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale/collator.hpp>
#include <boost/locale/generator.hpp>
#include <boost/locale/localization_backend.hpp>
#include "boostLocale/test/tools.hpp"
#include "boostLocale/test/unit_test.hpp"
#include <iomanip>

template<typename Char>
void test_comp(const std::locale& l, // LCOV_EXCL_LINE
               const std::basic_string<Char>& left,
               const std::basic_string<Char>& right,
               const boost::locale::collate_level level,
               const int expected)
{
    typedef std::basic_string<Char> string_type;
    TEST_EQ(boost::locale::comparator<Char>(l, level)(left, right), (expected < 0));
    if(level == boost::locale::collate_level::identical) {
        const std::collate<Char>& coll = std::use_facet<std::collate<Char>>(l);
        string_type lt = coll.transform(left.c_str(), left.c_str() + left.size());
        string_type rt = coll.transform(right.c_str(), right.c_str() + right.size());
        if(expected < 0)
            TEST_LT(lt, rt);
        else if(expected == 0)
            TEST_EQ(lt, rt);
        else
            TEST_GT(lt, rt);
        long lh = coll.hash(left.c_str(), left.c_str() + left.size());
        long rh = coll.hash(right.c_str(), right.c_str() + right.size());
        if(expected == 0)
            TEST_EQ(lh, rh);
        else
            TEST_NE(lh, rh);
    }
    const boost::locale::collator<Char>& coll = std::use_facet<boost::locale::collator<Char>>(l);
    string_type lt = coll.transform(level, left.c_str(), left.c_str() + left.size());
    TEST_EQ(lt, coll.transform(level, left));
    string_type rt = coll.transform(level, right.c_str(), right.c_str() + right.size());
    TEST_EQ(rt, coll.transform(level, right));
    if(expected < 0)
        TEST_LT(lt, rt);
    else if(expected == 0)
        TEST_EQ(lt, rt);
    else
        TEST_GT(lt, rt);
    long lh = coll.hash(level, left.c_str(), left.c_str() + left.size());
    TEST_EQ(lh, coll.hash(level, left));
    long rh = coll.hash(level, right.c_str(), right.c_str() + right.size());
    TEST_EQ(rh, coll.hash(level, right));
    if(expected == 0)
        TEST_EQ(lh, rh);
    else
        TEST_NE(lh, rh);
} // LCOV_EXCL_LINE

#define TEST_COMP(c, _l, _r) test_comp<c>(l, _l, _r, level, expected)

void compare(const std::string left, // LCOV_EXCL_LINE
             const std::string right,
             const boost::locale::collate_level level,
             const int expected)
{
    boost::locale::generator gen;
    std::locale l = gen("en_US.UTF-8");
    if(level == boost::locale::collate_level::identical)
        TEST_EQ(l(left, right), (expected < 0));
    TEST_COMP(char, left, right);
    TEST_COMP(wchar_t, to<wchar_t>(left), to<wchar_t>(right));
#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
    TEST_COMP(char16_t, to<char16_t>(left), to<char16_t>(right));
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
    TEST_COMP(char32_t, to<char32_t>(left), to<char32_t>(right));
#endif
} // LCOV_EXCL_LINE

void test_collate() // LCOV_EXCL_LINE
{
    constexpr int le = -1, gt = 1, eq = 0;
    using boost::locale::collate_level;

    compare("a", "A", collate_level::primary, eq);
    compare("a", "A", collate_level::secondary, eq);
    compare("A", "a", collate_level::tertiary, gt);
    compare("a", "A", collate_level::tertiary, le);
    compare("a", "A", collate_level::quaternary, le);
    compare("A", "a", collate_level::quaternary, gt);
    compare("a", "A", collate_level::identical, le);
    compare("A", "a", collate_level::identical, gt);
    compare("a", "ä", collate_level::primary, eq);
    compare("a", "ä", collate_level::secondary, le);
    compare("ä", "a", collate_level::secondary, gt);
    compare("a", "ä", collate_level::quaternary, le);
    compare("ä", "a", collate_level::quaternary, gt);
    compare("a", "ä", collate_level::identical, le);
    compare("ä", "a", collate_level::identical, gt);
    compare("a", "a", collate_level::identical, eq);
    compare("ä", "ä", collate_level::identical, eq);
} // LCOV_EXCL_LINE

BOOST_LOCALE_DISABLE_UNREACHABLE_CODE_WARNING
void test_main(int /*argc*/, char** /*argv*/)
{
#ifdef BOOST_LOCALE_NO_WINAPI_BACKEND
    std::cout << "WinAPI Backend is not build... Skipping\n";
    return;
#endif
    boost::locale::localization_backend_manager mgr = boost::locale::localization_backend_manager::global();
    mgr.select("winapi");
    boost::locale::localization_backend_manager::global(mgr);

    test_collate();
}

// boostinspect:noascii
