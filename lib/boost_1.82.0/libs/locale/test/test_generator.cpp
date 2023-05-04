//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale.hpp>
#include "../src/boost/locale/win32/lcid.hpp"
#include "boostLocale/test/tools.hpp"
#include "boostLocale/test/unit_test.hpp"
#include <boost/core/ignore_unused.hpp>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <vector>
#ifdef BOOST_LOCALE_WITH_ICU
#    include <unicode/uversion.h>
#    define BOOST_LOCALE_ICU_VERSION (U_ICU_VERSION_MAJOR_NUM * 100 + U_ICU_VERSION_MINOR_NUM)
#else
#    define BOOST_LOCALE_ICU_VERSION 0
#endif

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

namespace bl = boost::locale;

bool hasLocaleForBackend(const std::string& locale_name, const std::string& backendName)
{
    if(backendName == "winapi") {
#ifdef BOOST_LOCALE_NO_WINAPI_BACKEND
        boost::ignore_unused(locale_name);
        return false;
#else
        return bl::impl_win::locale_to_lcid(locale_name) != 0;
#endif
    } else if(backendName == "std")
        return has_std_locale(locale_name.c_str());
    else if(backendName == "posix")
        return has_posix_locale(locale_name);
    else                                         // ICU
        return BOOST_LOCALE_ICU_VERSION >= 5901; // First version to use (correct) CLDR data
}

void test_special_locales(const std::vector<std::string>& backends)
{
    for(const std::string& backendName : backends) {
        std::cout << "Backend: " << backendName << std::endl;
        bl::localization_backend_manager tmp_backend = bl::localization_backend_manager::global();
        tmp_backend.select(backendName);
        bl::localization_backend_manager::global(tmp_backend);

        bl::generator g;
        namespace as = bl::as;
        constexpr time_t datetime = 60 * 60 * 24 * (31 + 4) // Feb 5th
                                    + (15 * 60 + 42) * 60;  // 15:42

        const std::string enWorldName = "en_001.UTF-8";
        if(!hasLocaleForBackend(enWorldName, backendName))
            std::cout << "\tSkipping due to missing locale " << enWorldName << std::endl;
        else {
            auto l = g(enWorldName);
            const auto& info = std::use_facet<bl::info>(l);
            TEST_EQ(info.language(), "en");
            TEST_EQ(info.country(), "001");
            TEST(info.utf8());
            TEST_EQ(info.encoding(), "UTF-8");

            std::ostringstream os;
            os.imbue(l);
            os << as::time << as::gmt << as::time_short;
            os << datetime;
            TEST_EQ(os.str().substr(0, 4), "3:42"); // 3:42 pm
        }
        const std::string enEuropeName = "en_150.UTF-8";
        if(!hasLocaleForBackend(enEuropeName, backendName))
            std::cout << "\tSkipping due to missing locale " << enEuropeName << std::endl;
        else {
            auto l = g(enEuropeName);
            const auto& info = std::use_facet<bl::info>(l);
            TEST_EQ(info.language(), "en");
            TEST_EQ(info.country(), "150");
            TEST(info.utf8());
            TEST_EQ(info.encoding(), "UTF-8");

            std::ostringstream os;

            std::string expectedTimeFormat = "15:42";
            // The std locale may not fully support the 150 region and use a different format
            if(backendName == "std") {
                os.imbue(std::locale(os.getloc(), new std::time_put_byname<char>(enEuropeName)));
                os << std::put_time(gmtime_wrap(&datetime), "%X");
                expectedTimeFormat = os.str();
                os.str("");
            }

            os.imbue(l);
            os << as::time << as::gmt << as::time_short;
            os << datetime;
            TEST_EQ(os.str().substr(0, expectedTimeFormat.size()), expectedTimeFormat);
        }
    }
}

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

    const bl::localization_backend_manager orig_backend = bl::localization_backend_manager::global();
    for(const std::string& backendName : backends) {
        std::cout << "Backend: " << backendName << std::endl;
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
    std::cout << "Test special locales" << std::endl;
    test_special_locales(backends);
    std::cout << "Restoring original backend" << std::endl;
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
    TEST_EQ(std::use_facet<bl::info>(l).language(), "en");
    TEST_EQ(std::use_facet<bl::info>(l).country(), "US");
    TEST(!std::use_facet<bl::info>(l).utf8());
    TEST_EQ(std::use_facet<bl::info>(l).encoding(), "ISO8859-1");

    l = g("en_US.UTF-8");
    TEST_EQ(std::use_facet<bl::info>(l).language(), "en");
    TEST_EQ(std::use_facet<bl::info>(l).country(), "US");
    TEST(std::use_facet<bl::info>(l).utf8());
    TEST_EQ(std::use_facet<bl::info>(l).encoding(), "UTF-8");

    l = g("en_US.ISO8859-1");
    TEST_EQ(std::use_facet<bl::info>(l).language(), "en");
    TEST_EQ(std::use_facet<bl::info>(l).country(), "US");
    TEST(!std::use_facet<bl::info>(l).utf8());
    TEST_EQ(std::use_facet<bl::info>(l).encoding(), "ISO8859-1");

    // Check that generate() extends the given locale, not replaces it
    std::locale l_wt(std::locale::classic(), new test_facet);
    TEST(std::has_facet<test_facet>(g.generate(l_wt, "en_US.UTF-8")));
    TEST(!std::has_facet<test_facet>(g.generate("en_US.UTF-8")));
    TEST(std::has_facet<test_facet>(g.generate(l_wt, "en_US.ISO8859-1")));
    TEST(!std::has_facet<test_facet>(g.generate("en_US.ISO8859-1")));

    // Check caching works
    g.locale_cache_enabled(true);
    // Generate a locale with a specific facet which is then cached
    g.generate(l_wt, "en_US.UTF-8");
    g.generate(l_wt, "en_US.ISO8859-1");
    // Cached locale is returned -> facet is still there
    TEST(std::has_facet<test_facet>(g("en_US.UTF-8")));
    TEST(std::has_facet<test_facet>(g("en_US.ISO8859-1")));
    // Check a property to verify it doesn't simply return the same locale for each call
    TEST(std::use_facet<bl::info>(g("en_US.UTF-8")).utf8());
    TEST(!std::use_facet<bl::info>(g("en_US.ISO8859-1")).utf8());
}
