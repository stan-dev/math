//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale/encoding.hpp>
#include <boost/locale/formatting.hpp>
#include <boost/locale/generator.hpp>
#include <boost/locale/info.hpp>
#include <boost/locale/localization_backend.hpp>
#include <boost/core/ignore_unused.hpp>
#include <ctime>
#include <iomanip>
#include <iostream>
#ifndef BOOST_LOCALE_NO_POSIX_BACKEND
#    include <langinfo.h>
#    include <monetary.h>
#endif
#include "boostLocale/test/tools.hpp"
#include "boostLocale/test/unit_test.hpp"

#ifdef BOOST_LOCALE_NO_POSIX_BACKEND
// Dummy just to make it compile
size_t strftime_l(char*, size_t, const char*, const std::tm*, locale_t)
{
    return 0; // LCOV_EXCL_LINE
}
#endif

template<typename CharType>
std::basic_string<CharType> from_narrow(const std::string& s, locale_t lc)
{
#ifdef BOOST_LOCALE_NO_POSIX_BACKEND
    const std::string charset; // LCOV_EXCL_LINE
    boost::ignore_unused(lc);  // LCOV_EXCL_LINE
#else
    const std::string charset = nl_langinfo_l(CODESET, lc);
#endif
    return boost::locale::conv::to_utf<CharType>(s, charset);
}

template<>
std::basic_string<char> from_narrow(const std::string& s, locale_t)
{
    return s;
}

template<typename CharType>
void test_by_char(const std::locale& l, locale_t lreal)
{
    typedef std::basic_stringstream<CharType> ss_type;

    using namespace boost::locale;

    {
        std::cout << "- Testing as::posix" << std::endl;
        ss_type ss;
        ss.imbue(l);

        TEST(ss << 1045.45);
        double n;
        TEST(ss >> n);
        TEST_EQ(n, 1045.45);
        TEST_EQ(ss.str(), ascii_to<CharType>("1045.45"));
    }

    {
        std::cout << "- Testing as::number" << std::endl;
        ss_type ss;
        ss.imbue(l);

        ss << as::number;
        TEST(ss << 1045.45);
        double n;
        TEST(ss >> n);
        TEST_EQ(n, 1045.45);

        if(std::use_facet<boost::locale::info>(l).country() == "US")
            TEST_EQ(ss.str(), from_narrow<CharType>("1,045.45", lreal));
    }

    {
        std::cout << "- Testing as::currency national " << std::endl;

        char buf[64]{};
#ifndef BOOST_LOCALE_NO_POSIX_BACKEND
        TEST_GT(strfmon_l(buf, sizeof(buf), lreal, "%n", 1043.34), 0);
#endif

        ss_type ss;
        ss.imbue(l);

        ss << as::currency;
        TEST(ss << 1043.34);

        TEST_EQ(ss.str(), from_narrow<CharType>(buf, lreal));
    }

    {
        std::cout << "- Testing as::currency iso" << std::endl;
        char buf[64]{};
#ifndef BOOST_LOCALE_NO_POSIX_BACKEND
        TEST_GT(strfmon_l(buf, sizeof(buf), lreal, "%i", 1043.34), 0);
#endif
        ss_type ss;
        ss.imbue(l);

        ss << as::currency << as::currency_iso;
        TEST(ss << 1043.34);

        TEST_EQ(ss.str(), from_narrow<CharType>(buf, lreal));
    }

    {
        std::cout << "- Testing as::date/time" << std::endl;

        const time_t a_date = 3600 * 24 * (31 + 4);     // Feb 5th
        const time_t a_time = 3600 * 15 + 60 * 33 + 13; // 15:33:13
        const time_t a_datetime = a_date + a_time;

        char buf[64]{};
        const std::tm tm = *gmtime_wrap(&a_datetime);
        TEST_REQUIRE(strftime_l(buf, sizeof(buf), "%x", &tm, lreal) != 0u);
        const auto expDate = from_narrow<CharType>(buf, lreal);
        TEST_REQUIRE(strftime_l(buf, sizeof(buf), "%X", &tm, lreal) != 0u);
        const auto expTime = from_narrow<CharType>(buf, lreal);
        TEST_REQUIRE(strftime_l(buf, sizeof(buf), "%c", &tm, lreal) != 0u);
        const auto expDateTime = from_narrow<CharType>(buf, lreal);

        ss_type ss;
        ss.imbue(l);
        ss << as::time_zone("GMT");

        empty_stream(ss) << as::date << a_datetime;
        TEST_EQ(ss.str(), expDate);
        empty_stream(ss) << as::time << a_datetime;
        TEST_EQ(ss.str(), expTime);
        empty_stream(ss) << as::datetime << a_datetime;
        TEST_EQ(ss.str(), expDateTime);
        empty_stream(ss) << as::time_zone("GMT+01:00") << as::ftime(ascii_to<CharType>("%H")) << a_datetime;
        TEST_EQ(ss.str(), ascii_to<CharType>("16"));
        empty_stream(ss) << as::time_zone("GMT+00:15") << as::ftime(ascii_to<CharType>("%M")) << a_datetime;
        TEST_EQ(ss.str(), ascii_to<CharType>("48"));
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
    boost::locale::generator gen;
    for(const std::string locale_name : {"en_US.UTF-8", "en_US.ISO8859-1", "he_IL.UTF-8", "he_IL.ISO8859-8"}) {
        std::cout << locale_name << " locale" << std::endl;
        if(!has_posix_locale(locale_name))
            std::cout << locale_name << " not supported, skipping" << std::endl; // LCOV_EXCL_LINE
        else {
            std::locale generated_locale = gen(locale_name);
            locale_holder real_locale(newlocale(LC_ALL_MASK, locale_name.c_str(), nullptr));
            TEST_REQUIRE(real_locale);

            std::cout << "UTF-8" << std::endl;
            test_by_char<char>(generated_locale, real_locale);

            std::cout << "Wide UTF-" << sizeof(wchar_t) * 8 << std::endl;
            test_by_char<wchar_t>(generated_locale, real_locale);
        }
    }
    {
        std::cout << "Testing UTF-8 punct issues" << std::endl;
        const std::string locale_name = "ru_RU.UTF-8";
        if(!has_posix_locale(locale_name))
            std::cout << "- No Russian locale, skipping" << std::endl; // LCOV_EXCL_LINE
        else {
            std::ostringstream ss;
            ss.imbue(gen(locale_name));
            ss << std::setprecision(10) << boost::locale::as::number << 12345.45;
            const std::string v = ss.str();
            TEST(v == "12345,45" || v == "12 345,45" || v == "12.345,45");
        }
    }
}

// boostinspect:noascii
