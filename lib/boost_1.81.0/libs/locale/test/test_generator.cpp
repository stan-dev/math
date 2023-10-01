//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale.hpp>
#include "boostLocale/test/unit_test.hpp"
#include <iomanip>
#include <string>
#include <vector>

bool has_message(const std::locale& l)
{
    return std::has_facet<boost::locale::message_format<char>>(l);
}

struct test_facet : public std::locale::facet {
    test_facet() : std::locale::facet(0) {}
    static std::locale::id id;
};

std::locale::id test_facet::id;

template<typename CharType>
using codecvt_by_char_type = std::codecvt<CharType, char, std::mbstate_t>;

void test_main(int /*argc*/, char** /*argv*/)
{
    std::vector<std::string> backends;
#ifdef BOOST_LOCALE_WITH_ICU
    backends.push_back("icu");
#endif
#ifndef BOOST_LOCALE_NO_STD_BACKEND
    backends.push_back("std");
#endif
#ifndef BOOST_LOCALE_NO_WINAPI_BACKEND
    backends.push_back("winapi");
#endif
#ifndef BOOST_LOCALE_NO_POSIX_BACKEND
    backends.push_back("posix");
#endif

    namespace bl = boost::locale;
    const bl::localization_backend_manager orig_backend = bl::localization_backend_manager::global();
    for(const std::string& backendName : backends) {
        std::cout << "Backend: " << backendName << '\n';
        bl::localization_backend_manager tmp_backend = bl::localization_backend_manager::global();
        tmp_backend.select(backendName);
        bl::localization_backend_manager::global(tmp_backend);
        bl::generator g;
        const std::locale l = g("en_US.UTF-8");
#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
#    define TEST_HAS_FACET_CHAR16(facet, l) TEST(std::has_facet<facet<char16_t>>(l))
#else
#    define TEST_HAS_FACET_CHAR16(facet, l) (void)0
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
#    define TEST_HAS_FACET_CHAR32(facet, l) TEST(std::has_facet<facet<char32_t>>(l))
#else
#    define TEST_HAS_FACET_CHAR32(facet, l) (void)0
#endif
#define TEST_HAS_FACETS(facet, l)            \
    TEST(std::has_facet<facet<char>>(l));    \
    TEST(std::has_facet<facet<wchar_t>>(l)); \
    TEST_HAS_FACET_CHAR16(facet, l);         \
    TEST_HAS_FACET_CHAR32(facet, l)

        // Convert
        TEST_HAS_FACETS(bl::converter, l);
        TEST_HAS_FACETS(std::collate, l);
        // Formatting
        TEST_HAS_FACETS(std::num_put, l);
        TEST_HAS_FACETS(std::time_put, l);
        TEST_HAS_FACETS(std::numpunct, l);
        TEST_HAS_FACETS(std::moneypunct, l);
        // Parsing
        TEST_HAS_FACETS(std::num_get, l);
        // Message
        TEST_HAS_FACETS(bl::message_format, l);
        // Codepage
        TEST_HAS_FACETS(codecvt_by_char_type, l);
        // Boundary
        if(backendName == "icu") {
            TEST_HAS_FACETS(bl::boundary::boundary_indexing, l);
        }
        // calendar
        TEST(std::has_facet<bl::calendar_facet>(l));
        // information
        TEST(std::has_facet<bl::info>(l));
    }
    bl::localization_backend_manager::global(orig_backend);

    bl::generator g;
    std::locale l = g("en_US.UTF-8");
    TEST(has_message(l));
    g.categories(g.categories() ^ bl::category_t::message);
    g.locale_cache_enabled(true);
    g("en_US.UTF-8");
    g.categories(g.categories() | bl::category_t::message);
    l = g("en_US.UTF-8");
    TEST(!has_message(l));
    g.clear_cache();
    g.locale_cache_enabled(false);
    l = g("en_US.UTF-8");
    TEST(has_message(l));
    g.characters(g.characters() ^ bl::char_facet_t::char_f);
    l = g("en_US.UTF-8");
    TEST(!has_message(l));
    g.characters(g.characters() | bl::char_facet_t::char_f);
    l = g("en_US.UTF-8");
    TEST(has_message(l));

    l = g("en_US.ISO8859-1");
    TEST(std::use_facet<bl::info>(l).language() == "en");
    TEST(std::use_facet<bl::info>(l).country() == "US");
    TEST(!std::use_facet<bl::info>(l).utf8());
    TEST(std::use_facet<bl::info>(l).encoding() == "iso8859-1");

    l = g("en_US.UTF-8");
    TEST(std::use_facet<bl::info>(l).language() == "en");
    TEST(std::use_facet<bl::info>(l).country() == "US");
    TEST(std::use_facet<bl::info>(l).utf8());

    l = g("en_US.ISO8859-1");
    TEST(std::use_facet<bl::info>(l).language() == "en");
    TEST(std::use_facet<bl::info>(l).country() == "US");
    TEST(!std::use_facet<bl::info>(l).utf8());
    TEST(std::use_facet<bl::info>(l).encoding() == "iso8859-1");

    l = g("en_US.ISO8859-1");
    TEST(std::use_facet<bl::info>(l).language() == "en");
    TEST(std::use_facet<bl::info>(l).country() == "US");
    TEST(!std::use_facet<bl::info>(l).utf8());
    TEST(std::use_facet<bl::info>(l).encoding() == "iso8859-1");

    std::locale l_wt(std::locale::classic(), new test_facet);

    TEST(std::has_facet<test_facet>(g.generate(l_wt, "en_US.UTF-8")));
    TEST(std::has_facet<test_facet>(g.generate(l_wt, "en_US.ISO8859-1")));
    TEST(!std::has_facet<test_facet>(g("en_US.UTF-8")));
    TEST(!std::has_facet<test_facet>(g("en_US.ISO8859-1")));

    g.locale_cache_enabled(true);
    g.generate(l_wt, "en_US.UTF-8");
    g.generate(l_wt, "en_US.ISO8859-1");
    TEST(std::has_facet<test_facet>(g("en_US.UTF-8")));
    TEST(std::has_facet<test_facet>(g("en_US.ISO8859-1")));
    TEST(std::use_facet<bl::info>(g("en_US.UTF-8")).utf8());
    TEST(!std::use_facet<bl::info>(g("en_US.ISO8859-1")).utf8());
}

// boostinspect:noascii
