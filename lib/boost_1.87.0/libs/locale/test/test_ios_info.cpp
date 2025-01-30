//
// Copyright (c) 2022 Alexander Grund
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale/formatting.hpp>
#include "boostLocale/test/tools.hpp"
#include "boostLocale/test/unit_test.hpp"
#include <locale>
#include <sstream>
#include <stdexcept>

void test_member_methods()
{
    std::stringstream ss;
    {
        // Test the default constructed instance
        boost::locale::time_zone::global("Global TZ");
        const auto& info = boost::locale::ios_info::get(ss);
        TEST_EQ(info.display_flags(), 0u);
        TEST_EQ(info.currency_flags(), 0u);
        TEST_EQ(info.date_flags(), 0u);
        TEST_EQ(info.time_flags(), 0u);
        TEST_EQ(info.domain_id(), 0);
        TEST_EQ(info.time_zone(), "Global TZ");
        TEST_THROWS(info.date_time_pattern<char>(), std::bad_cast);
    }
    {
        namespace flags = boost::locale::flags;
        auto& info = boost::locale::ios_info::get(ss);
        for(const auto flag : {flags::posix,
                               flags::number,
                               flags::currency,
                               flags::percent,
                               flags::date,
                               flags::time,
                               flags::datetime,
                               flags::strftime,
                               flags::spellout,
                               flags::ordinal})
        {
            info.display_flags(flag);
            TEST_EQ(info.display_flags(), flag);
        }
        for(const auto flag : {flags::currency_default, flags::currency_iso, flags::currency_national}) {
            info.currency_flags(flag);
            TEST_EQ(info.currency_flags(), flag);
        }
        for(const auto flag :
            {flags::time_default, flags::time_short, flags::time_medium, flags::time_long, flags::time_full}) {
            info.time_flags(flag);
            TEST_EQ(info.time_flags(), flag);
        }
        for(const auto flag :
            {flags::date_default, flags::date_short, flags::date_medium, flags::date_long, flags::date_full}) {
            info.date_flags(flag);
            TEST_EQ(info.date_flags(), flag);
        }
        // Should only change 1 setting
        info.display_flags(flags::number);
        info.currency_flags(flags::currency_iso);
        info.time_flags(flags::time_medium);
        info.date_flags(flags::date_short);
        TEST_EQ(info.display_flags(), flags::number);
        TEST_EQ(info.currency_flags(), flags::currency_iso);
        TEST_EQ(info.time_flags(), flags::time_medium);
        TEST_EQ(info.date_flags(), flags::date_short);

        info.display_flags(flags::ordinal);
        TEST_EQ(info.display_flags(), flags::ordinal);
        TEST_EQ(info.currency_flags(), flags::currency_iso);
        TEST_EQ(info.time_flags(), flags::time_medium);
        TEST_EQ(info.date_flags(), flags::date_short);

        info.display_flags(flags::number);
        info.currency_flags(flags::currency_national);
        TEST_EQ(info.display_flags(), flags::number);
        TEST_EQ(info.currency_flags(), flags::currency_national);
        TEST_EQ(info.time_flags(), flags::time_medium);
        TEST_EQ(info.date_flags(), flags::date_short);

        info.display_flags(flags::number);
        info.currency_flags(flags::currency_iso);
        info.time_flags(flags::time_full);
        TEST_EQ(info.display_flags(), flags::number);
        TEST_EQ(info.currency_flags(), flags::currency_iso);
        TEST_EQ(info.time_flags(), flags::time_full);
        TEST_EQ(info.date_flags(), flags::date_short);

        info.display_flags(flags::number);
        info.currency_flags(flags::currency_iso);
        info.time_flags(flags::time_medium);
        info.date_flags(flags::date_full);
        TEST_EQ(info.display_flags(), flags::number);
        TEST_EQ(info.currency_flags(), flags::currency_iso);
        TEST_EQ(info.time_flags(), flags::time_medium);
        TEST_EQ(info.date_flags(), flags::date_full);

        info.domain_id(42);
        TEST_EQ(info.domain_id(), 42);

        info.time_zone("Test-TZ");
        TEST_EQ(info.time_zone(), "Test-TZ");

        info.date_time_pattern(std::string("Pattern"));
        TEST_EQ(info.date_time_pattern<char>(), "Pattern");
    }
}

void test_manipulators()
{
    std::stringstream ss;
    const auto& info = boost::locale::ios_info::get(ss);
    namespace as = boost::locale::as;
    namespace flags = boost::locale::flags;
    ss << as::number << as::currency_iso << as::time_short << as::date_short << as::time_zone("FooTZ");
    TEST_EQ(info.display_flags(), flags::number);
    TEST_EQ(info.currency_flags(), flags::currency_iso);
    TEST_EQ(info.time_flags(), flags::time_short);
    TEST_EQ(info.date_flags(), flags::date_short);
    TEST_EQ(info.time_zone(), "FooTZ");
    TEST_THROWS(info.date_time_pattern<char>(), std::bad_cast);

#define TEST_MODIFIER(modifier) \
    ss << as::modifier;         \
    TEST_EQ(info.display_flags(), flags::modifier);

    TEST_MODIFIER(posix);
    TEST_MODIFIER(currency);
    TEST_MODIFIER(percent);
    TEST_MODIFIER(date);
    TEST_MODIFIER(time);
    TEST_MODIFIER(strftime);
    TEST_MODIFIER(spellout);
    TEST_MODIFIER(ordinal);
#undef TEST_MODIFIER

#define TEST_MODIFIER(modifier) \
    ss << as::modifier;         \
    TEST_EQ(info.currency_flags(), flags::modifier);

    TEST_MODIFIER(currency_default);
    TEST_MODIFIER(currency_iso);
    TEST_MODIFIER(currency_national);
#undef TEST_MODIFIER

#define TEST_MODIFIER(modifier) \
    ss << as::modifier;         \
    TEST_EQ(info.currency_flags(), flags::modifier);

    TEST_MODIFIER(currency_default);
    TEST_MODIFIER(currency_iso);
    TEST_MODIFIER(currency_national);
#undef TEST_MODIFIER

#define TEST_MODIFIER(modifier) \
    ss << as::modifier;         \
    TEST_EQ(info.time_flags(), flags::modifier);

    TEST_MODIFIER(time_default);
    TEST_MODIFIER(time_short);
    TEST_MODIFIER(time_medium);
    TEST_MODIFIER(time_long);
    TEST_MODIFIER(time_full);
#undef TEST_MODIFIER

#define TEST_MODIFIER(modifier) \
    ss << as::modifier;         \
    TEST_EQ(info.date_flags(), flags::modifier);

    TEST_MODIFIER(date_default);
    TEST_MODIFIER(date_short);
    TEST_MODIFIER(date_medium);
    TEST_MODIFIER(date_long);
    TEST_MODIFIER(date_full);
#undef TEST_MODIFIER

    ss << as::number << as::ftime("Test format");
    TEST_EQ(info.display_flags(), flags::strftime); // Auto changed
    TEST_EQ(info.date_time_pattern<char>(), "Test format");
    ss >> as::number >> as::ftime("Test format2");  // Also as input
    TEST_EQ(info.display_flags(), flags::strftime); // Auto changed
    TEST_EQ(info.date_time_pattern<char>(), "Test format2");

    ss << as::gmt;
    TEST_EQ(info.time_zone(), "GMT");
    boost::locale::time_zone::global("Global TZ");
    ss << as::local_time;
    TEST_EQ(info.time_zone(), "Global TZ");
    ss << as::time_zone("Foo TZ");
    TEST_EQ(info.time_zone(), "Foo TZ");
    ss >> as::time_zone("Bar TZ"); // Also as input
    TEST_EQ(info.time_zone(), "Bar TZ");

    std::wistringstream ss2;
    ss2 >> as::ftime(L"My TZ");
    const auto& info2 = boost::locale::ios_info::get(ss2);
    TEST_EQ(info2.display_flags(), flags::strftime);
    TEST_EQ(info2.date_time_pattern<wchar_t>(), L"My TZ");
}

void test_any_string()
{
    boost::locale::detail::any_string s;
    TEST_THROWS(s.get<char>(), std::bad_cast);
    TEST_THROWS(s.get<wchar_t>(), std::bad_cast);

    s.set<char>("Char Pattern");
    TEST_EQ(s.get<char>(), "Char Pattern");
    TEST_THROWS(s.get<wchar_t>(), std::bad_cast);

    s.set<wchar_t>(ascii_to<wchar_t>("WChar Pattern"));
    TEST_EQ(s.get<wchar_t>(), ascii_to<wchar_t>("WChar Pattern"));
    TEST_THROWS(s.get<char>(), std::bad_cast);

    s.set<char16_t>(ascii_to<char16_t>("Char16 Pattern"));
    TEST_EQ(s.get<char16_t>(), ascii_to<char16_t>("Char16 Pattern"));
    TEST_THROWS(s.get<char>(), std::bad_cast);

    s.set<char32_t>(ascii_to<char32_t>("Char32 Pattern"));
    TEST_EQ(s.get<char32_t>(), ascii_to<char32_t>("Char32 Pattern"));
    TEST_THROWS(s.get<char16_t>(), std::bad_cast);

#ifndef BOOST_LOCALE_NO_CXX20_STRING8
    s.set<char8_t>(ascii_to<char8_t>("Char8 Pattern"));
    TEST_EQ(s.get<char8_t>(), ascii_to<char8_t>("Char8 Pattern"));
    TEST_THROWS(s.get<char32_t>(), std::bad_cast);
#endif

    boost::locale::detail::any_string s1, s2, empty;
    s1.set<char>("Char");
    s2.set<wchar_t>(ascii_to<wchar_t>("WChar"));
    // Copy ctor
    boost::locale::detail::any_string s3(s1);
    TEST_EQ(s3.get<char>(), "Char");
    TEST_EQ(s1.get<char>(), "Char");
    // Ensure deep copy
    s3.set<char>("Foo");
    TEST_EQ(s3.get<char>(), "Foo");
    TEST_EQ(s1.get<char>(), "Char");
    // Copy assign
    s3 = s2;
    TEST_EQ(s3.get<wchar_t>(), ascii_to<wchar_t>("WChar"));
    TEST_EQ(s2.get<wchar_t>(), ascii_to<wchar_t>("WChar"));
    // Move assign
    s3 = std::move(s1);
    TEST_EQ(s3.get<char>(), "Char");
    // From empty
    s3 = empty;
    TEST_THROWS(s3.get<char>(), std::bad_cast);
}

void test_main(int /*argc*/, char** /*argv*/)
{
    test_any_string();
    test_member_methods();
    test_manipulators();
}
