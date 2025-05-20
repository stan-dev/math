//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
// Copyright (c) 2021-2023 Alexander Grund
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale.hpp>
#include <boost/locale/conversion.hpp>
#include "../src/boost/locale/win32/lcid.hpp"
#include "boostLocale/test/tools.hpp"
#include "boostLocale/test/unit_test.hpp"
#include <boost/assert.hpp>
#include <boost/core/ignore_unused.hpp>
#include <algorithm>
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

namespace boost { namespace locale { namespace test {
    template<class Facet>
    BOOST_NOINLINE bool is_facet(const std::locale::facet* facet)
    {
        return dynamic_cast<const Facet*>(facet) != nullptr;
    }

    template<class Facet>
    bool has_facet(const std::locale& l)
    {
        return std::has_facet<Facet>(l) && is_facet<Facet>(&std::use_facet<Facet>(l));
    }

    template<class Facet>
    bool has_not_facet(const std::locale& l)
    {
        const Facet* f;
        try {
            f = &std::use_facet<Facet>(l);
        } catch(const std::bad_cast&) {
            return !std::has_facet<Facet>(l);
        }
        // This mustn't be reached, checks for debugging
        TEST(is_facet<Facet>(f));        // LCOV_EXCL_LINE
        TEST(!std::has_facet<Facet>(l)); // LCOV_EXCL_LINE
        return false;                    // LCOV_EXCL_LINE
    }
}}} // namespace boost::locale::test
namespace blt = boost::locale::test;

bool has_message(const std::locale& l)
{
    return blt::has_facet<boost::locale::message_format<char>>(l);
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
    if(backendName == "winapi")
        return has_win_locale(locale_name);
    else if(backendName == "std")
        return has_std_locale(locale_name.c_str());
    else if(backendName == "posix")
        return has_posix_locale(locale_name);
    else {
        BOOST_ASSERT(backendName == "icu");
        return BOOST_LOCALE_ICU_VERSION >= 5901; // First version to use (correct) CLDR data
    }
}

void test_special_locales()
{
    bl::localization_backend_manager backend = bl::localization_backend_manager::global();
    for(const std::string& backendName : backend.get_all_backends()) {
        std::cout << "Backend: " << backendName << std::endl;
        backend.select(backendName);
        bl::localization_backend_manager::global(backend);

        {
            const auto utf8LocaleName = bl::util::get_system_locale(true);
            // The WinAPI backend only supports UTF-8 encoding and hence always returns the UTF-8 locale
            const auto ansiLocaleName = (backendName == "winapi") ? utf8LocaleName : bl::util::get_system_locale(false);
            bl::generator g;
            g.use_ansi_encoding(true);
            std::locale l = g("");
            TEST_EQ(std::use_facet<bl::info>(l).name(), ansiLocaleName);
            g.use_ansi_encoding(false);
            l = g("");
            TEST_EQ(std::use_facet<bl::info>(l).name(), utf8LocaleName);
            g.use_ansi_encoding(true);
            l = g("");
            TEST_EQ(std::use_facet<bl::info>(l).name(), ansiLocaleName);
        }

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
                empty_stream(os) << std::put_time(gmtime_wrap(&datetime), "%X");
                expectedTimeFormat = os.str();
            }

            os.imbue(l);
            empty_stream(os) << as::time << as::gmt << as::time_short;
            os << datetime;
            TEST_EQ(os.str().substr(0, expectedTimeFormat.size()), expectedTimeFormat);
        }
    }
}

bool has_unicode_classic_locale()
{
    std::locale l = std::locale::classic();
    for(const auto name : {"C.UTF-8", "C.utf8"}) {
        try {
            l = std::locale(name);
            break;
        } catch(...) {
        }
    }
    const wchar_t s = L'\u03B4';
    // Check that that the Unicode character is handled
    return std::use_facet<std::ctype<wchar_t>>(l).toupper(s) != s;
}

// For a non-existing locale the C locale will be used as a fallback
// If UTF-8 is requested/reported then UTF-8 will still be used as much as possible
void test_invalid_locale()
{
    std::ostringstream classicStream;
    classicStream.imbue(std::locale::classic());
    const boost::locale::util::locale_data systemLocale(boost::locale::util::get_system_locale());

    bl::localization_backend_manager tmp_backend = bl::localization_backend_manager::global();
    tmp_backend.select("std");
    bl::localization_backend_manager::global(tmp_backend);
    bl::generator g;
    std::locale nonExistingLocale = g("noLang_noCountry." + systemLocale.encoding());
    const auto& info = std::use_facet<bl::info>(nonExistingLocale);
    TEST_EQ(info.language(), "nolang");
    TEST_EQ(info.country(), "NOCOUNTRY");
    std::ostringstream os;
    os.imbue(nonExistingLocale);
    os << boost::locale::as::number << 123456789 << " " << 1234567.89;
    classicStream << 123456789 << " " << 1234567.89;
    TEST_EQ(os.str(), classicStream.str());

    // Request UTF-8 explicitly and check that it is used even when the locale doesn't exist
    if(!info.utf8())
        nonExistingLocale = g("noLang_noCountry.UTF-8");
    if(has_unicode_classic_locale()) {
        // Case conversion works only if the backend supports classic locale with Unicode
        TEST_EQ(boost::locale::to_upper("δ", nonExistingLocale), "Δ");
    }
    // The codecvt facet always supports UTF-8
    {
        auto& cvt = std::use_facet<std::codecvt<wchar_t, char, std::mbstate_t>>(nonExistingLocale);
        // String with Unicode chars from different cultures
        const std::wstring wide_str = L"\U0001D49E-\u043F\u0440\u0438\u0432\u0435\u0442-\u3084\u3042";

        std::mbstate_t state{};
        const wchar_t* from_next = wide_str.c_str();
        const wchar_t* from_end = from_next + wide_str.size();
        char out_str[32]{};
        char* to_next = out_str;
        TEST_EQ(cvt.out(state, from_next, from_end, from_next, out_str, out_str + sizeof(out_str), to_next), cvt.ok);
        const std::string utf8_str = boost::locale::conv::utf_to_utf<char>(wide_str);
        TEST_EQ(out_str, utf8_str);
    }
}

void test_install_chartype(const std::string& backendName)
{
    // Use ASCII and UTF-8 encoding
    for(const std::string localeName : {"C", "en_US.UTF-8"}) {
        std::cout << "--- Locale: " << localeName << std::endl;
        const std::locale origLocale = bl::generator{}(localeName);
        const auto backend = bl::localization_backend_manager::global().create();
        backend->set_option("locale", localeName);
        for(auto category = bl::per_character_facet_first; category <= bl::per_character_facet_last; ++category) {
            std::cout << "---- Testing category " << static_cast<unsigned>(category) << '\n';
            // This should modify the locale
            const std::locale newLocale_char = backend->install(origLocale, category, bl::char_facet_t::char_f);
            // This should not
            const std::locale newLocale_nochar = backend->install(origLocale, category, bl::char_facet_t::nochar);
            // But the boundary facet is only implemented in ICU, so for all else the locale is still unchanged
            if(category != bl::category_t::boundary || backendName == "icu")
                TEST(origLocale != newLocale_char);
            else
                TEST(origLocale == newLocale_char);
            TEST(origLocale == newLocale_nochar);
        }
    }
}

template<typename Char>
struct dummy_collate : std::collate<Char> {};

template<typename Char>
bool has_dummy_collate(const std::locale& l)
{
    const auto& col = std::use_facet<std::collate<Char>>(l); // Implicitely require existance of std::collate
    return blt::is_facet<dummy_collate<Char>>(&col);
}

void test_std_collate_replaced(const std::string& /*backendName*/)
{
    std::locale origLocale = std::locale::classic();
    origLocale = std::locale(origLocale, new dummy_collate<char>);
    origLocale = std::locale(origLocale, new dummy_collate<wchar_t>);
#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
    origLocale = std::locale(origLocale, new dummy_collate<char16_t>);
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
    origLocale = std::locale(origLocale, new dummy_collate<char32_t>);
#endif

    // Use ASCII and UTF-8 encoding
    for(const std::string localeName : {"C", "en_US.UTF-8"}) {
        std::cout << "--- Locale: " << localeName << std::endl;
        bl::generator g;
        g.categories(boost::locale::category_t::collation);
        const std::locale l = g.generate(origLocale, localeName);
        TEST(has_dummy_collate<char>(origLocale));
        TEST(!has_dummy_collate<char>(l));
        TEST(has_dummy_collate<wchar_t>(origLocale));
        TEST(!has_dummy_collate<wchar_t>(l));
#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
        TEST(has_dummy_collate<char16_t>(origLocale));
        TEST(!has_dummy_collate<char16_t>(l));
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
        TEST(has_dummy_collate<char32_t>(origLocale));
        TEST(!has_dummy_collate<char32_t>(l));
#endif
    }
}

void test_main(int /*argc*/, char** /*argv*/)
{
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
        std::sort(backends.begin(), backends.end());

        std::vector<std::string> all_backends = bl::localization_backend_manager::global().get_all_backends();
        std::sort(all_backends.begin(), all_backends.end());
        TEST_EQ(all_backends, backends);
    }

    const bl::localization_backend_manager orig_backend = bl::localization_backend_manager::global();
    for(const std::string& backendName : orig_backend.get_all_backends()) {
        std::cout << "Backend: " << backendName << std::endl;
        bl::localization_backend_manager tmp_backend = bl::localization_backend_manager::global();
        tmp_backend.select(backendName);
        bl::localization_backend_manager::global(tmp_backend);
        bl::generator g;
        for(const std::string localeName : {"", "C", "en_US.UTF-8", "en_US.ISO8859-1", "tr_TR.windows1254"}) {
            std::cout << "-- Locale: " << localeName << std::endl;
            const std::locale l = g(localeName);
#ifdef __cpp_char8_t
#    define TEST_FOR_CHAR8(check) TEST(check)
#else
#    define TEST_FOR_CHAR8(check) (void)0
#endif
#ifndef BOOST_LOCALE_NO_CXX20_STRING8
#    define TEST_FOR_STRING8(check) TEST(check)
#else
#    define TEST_FOR_STRING8(check) (void)0
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
#    define TEST_FOR_CHAR16(check) TEST(check)
#else
#    define TEST_FOR_CHAR16(check) (void)0
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
#    define TEST_FOR_CHAR32(check) TEST(check)
#else
#    define TEST_FOR_CHAR32(check) (void)0
#endif
#define TEST_HAS_FACET_CHAR8(facet, l) TEST_FOR_CHAR8(blt::has_facet<facet<char8_t>>(l))
#define TEST_HAS_FACET_CHAR16(facet, l) TEST_FOR_CHAR16(blt::has_facet<facet<char16_t>>(l))
#define TEST_HAS_FACET_CHAR32(facet, l) TEST_FOR_CHAR32(blt::has_facet<facet<char32_t>>(l))
#define TEST_HAS_FACET_STRING8(facet, l) TEST_FOR_STRING8(blt::has_facet<facet<char8_t>>(l))

#define TEST_HAS_FACETS(facet, l)                \
    do {                                         \
        TEST(blt::has_facet<facet<char>>(l));    \
        TEST(blt::has_facet<facet<wchar_t>>(l)); \
        TEST_HAS_FACET_CHAR16(facet, l);         \
        TEST_HAS_FACET_CHAR32(facet, l);         \
    } while(false)

            // Convert
            TEST_HAS_FACETS(bl::converter, l);
            TEST_HAS_FACET_STRING8(bl::converter, l);
            // Collator
            TEST_HAS_FACETS(std::collate, l);
            if(backendName == "icu" || (backendName == "winapi" && std::use_facet<bl::info>(l).utf8())) {
                TEST_HAS_FACETS(bl::collator, l);
                TEST_HAS_FACET_STRING8(bl::collator, l);
            } else {
                TEST(blt::has_not_facet<bl::collator<char>>(l));
                TEST(blt::has_not_facet<bl::collator<wchar_t>>(l));
                TEST_FOR_STRING8(blt::has_not_facet<bl::collator<char8_t>>(l));
                TEST_FOR_CHAR16(blt::has_not_facet<bl::collator<char16_t>>(l));
                TEST_FOR_CHAR32(blt::has_not_facet<bl::collator<char32_t>>(l));
            }
            // Formatting
            TEST_HAS_FACETS(std::num_put, l);
            TEST_HAS_FACETS(std::time_put, l);
            TEST_HAS_FACETS(std::numpunct, l);
            TEST_HAS_FACETS(std::moneypunct, l);
            // Parsing
            TEST_HAS_FACETS(std::num_get, l);
            // Message
            TEST_HAS_FACETS(bl::message_format, l);
            TEST_HAS_FACET_STRING8(bl::message_format, l);
            // Codepage
            TEST_HAS_FACETS(codecvt_by_char_type, l);
            // Boundary
            if(backendName == "icu") {
                TEST_HAS_FACETS(bl::boundary::boundary_indexing, l);
                TEST_HAS_FACET_CHAR8(bl::boundary::boundary_indexing, l);
            }
            // calendar
            TEST(blt::has_facet<bl::calendar_facet>(l));
            // information
            TEST(blt::has_facet<bl::info>(l));
        }

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
        {
            const auto& info = std::use_facet<bl::info>(l);
            TEST_EQ(info.language(), "en");
            TEST_EQ(info.country(), "US");
            TEST(!info.utf8());
            TEST_EQ(info.variant(), "");
            TEST_EQ(info.encoding(), "ISO8859-1");
        }
        l = g("en_US.UTF-8");
        {
            const auto& info = std::use_facet<bl::info>(l);
            TEST_EQ(info.language(), "en");
            TEST_EQ(info.country(), "US");
            TEST(info.utf8());
            TEST_EQ(info.variant(), "");
            TEST_EQ(info.encoding(), "UTF-8");
        }
        l = g("da_DK.ISO8859-15@euro");
        {
            const auto& info = std::use_facet<bl::info>(l);
            TEST_EQ(info.language(), "da");
            TEST_EQ(info.country(), "DK");
            TEST(!info.utf8());
            TEST_EQ(info.variant(), "euro");
            TEST_EQ(info.encoding(), "ISO8859-15");
        }
        l = g("en_US.ISO8859-1");
        {
            const auto& info = std::use_facet<bl::info>(l);
            TEST_EQ(info.language(), "en");
            TEST_EQ(info.country(), "US");
            TEST(!info.utf8());
            TEST_EQ(info.variant(), "");
            TEST_EQ(info.encoding(), "ISO8859-1");
        }

        // Check that generate() extends the given locale, not replaces it
        std::locale l_wt(std::locale::classic(), new test_facet);
        TEST(blt::has_facet<test_facet>(g.generate(l_wt, "en_US.UTF-8")));
        TEST(!blt::has_facet<test_facet>(g.generate("en_US.UTF-8")));
        TEST(blt::has_facet<test_facet>(g.generate(l_wt, "en_US.ISO8859-1")));
        TEST(!blt::has_facet<test_facet>(g.generate("en_US.ISO8859-1")));

        // Check caching works
        g.locale_cache_enabled(true);
        // Generate a locale with a specific facet which is then cached
        g.generate(l_wt, "en_US.UTF-8");
        g.generate(l_wt, "en_US.ISO8859-1");
        // Cached locale is returned -> facet is still there
        TEST(blt::has_facet<test_facet>(g("en_US.UTF-8")));
        TEST(blt::has_facet<test_facet>(g("en_US.ISO8859-1")));
        // Check a property to verify it doesn't simply return the same locale for each call
        TEST(std::use_facet<bl::info>(g("en_US.UTF-8")).utf8());
        TEST(!std::use_facet<bl::info>(g("en_US.ISO8859-1")).utf8());

        test_install_chartype(backendName);
        test_std_collate_replaced(backendName);
    }
    std::cout << "Test special locales" << std::endl;
    test_special_locales();
    test_invalid_locale();
}
