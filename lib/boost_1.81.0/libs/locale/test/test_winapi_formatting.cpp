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
#ifndef BOOST_LOCALE_NO_WINAPI_BACKEND
#    ifndef NOMINMAX
#        define NOMINMAX
#    endif
#    include <windows.h>
#endif
#include "../src/boost/locale/win32/lcid.hpp"
#include "boostLocale/test/tools.hpp"
#include "boostLocale/test/unit_test.hpp"

template<typename CharType>
void test_by_char(const std::locale& l, std::string name, int lcid)
{
    typedef std::basic_stringstream<CharType> ss_type;

    using namespace boost::locale;

    {
        std::cout << "--- Testing as::posix" << std::endl;
        ss_type ss;
        ss.imbue(l);

        TEST(ss << 1045.45);
        double n;
        TEST(ss >> n);
        TEST(n == 1045.45);
        TEST_EQ(to_utf8(ss.str()), "1045.45");
    }

    {
        std::cout << "--- Testing as::number" << std::endl;
        ss_type ss;
        ss.imbue(l);

        ss << as::number;
        TEST(ss << 1045.45);
        double n;
        TEST(ss >> n);
        TEST_EQ(n, 1045.45);

        if(name == "ru_RU.UTF-8") {
            BOOST_LOCALE_START_CONST_CONDITION
            if(sizeof(CharType) == 1) // Space
                TEST_EQ(to_utf8(ss.str()), "1 045,45");
            else // Non-breaking space
                TEST_EQ(to_utf8(ss.str()),
                        "1\xC2\xA0"
                        "045,45");
            BOOST_LOCALE_END_CONST_CONDITION
        } else
            TEST_EQ(to_utf8(ss.str()), "1,045.45");
    }

    {
        std::cout << "--- Testing as::currency " << std::endl;

        ss_type ss;
        ss.imbue(l);

        ss << as::currency;
        TEST(ss << 1043.34);

#ifndef BOOST_LOCALE_NO_WINAPI_BACKEND
        wchar_t buf[256];
        GetCurrencyFormatW(lcid, 0, L"1043.34", 0, buf, 256);
        TEST_EQ(to_utf8(ss.str()), to_utf8<wchar_t>(buf));
#else
        boost::ignore_unused(lcid);
#endif
    }

    {
        std::cout << "--- Testing as::date/time" << std::endl;
        ss_type ss;
        ss.imbue(l);

        time_t a_date = 3600 * 24 * (31 + 4); // Feb 5th
        time_t a_time = 3600 * 15 + 60 * 33;  // 15:33:13
        time_t a_timesec = 13;
        time_t a_datetime = a_date + a_time + a_timesec;

        ss << as::time_zone("GMT");

        ss << as::date << a_datetime << CharType('\n');
        ss << as::time << a_datetime << CharType('\n');
        ss << as::datetime << a_datetime << CharType('\n');
        ss << as::time_zone("GMT+01:00");
        ss << as::ftime(ascii_to<CharType>("%H")) << a_datetime << CharType('\n');
        ss << as::time_zone("GMT+00:15");
        ss << as::ftime(ascii_to<CharType>("%M")) << a_datetime << CharType('\n');

#ifndef BOOST_LOCALE_NO_WINAPI_BACKEND
        wchar_t time_buf[256];
        wchar_t date_buf[256];
        SYSTEMTIME st = {1970, 2, 5, 5, 15, 33, 13, 0};
        GetTimeFormatW(lcid, 0, &st, 0, time_buf, 256);
        GetDateFormatW(lcid, 0, &st, 0, date_buf, 256);
        TEST_EQ(
          to_utf8(ss.str()),
          to_utf8(std::wstring(date_buf) + L"\n" + time_buf + L"\n" + date_buf + L" " + time_buf + L"\n16\n48\n"));
#endif
    }
}

void test_date_time(std::locale l)
{
    std::ostringstream ss;
    ss.imbue(l);

    ss << boost::locale::as::time_zone("GMT");

    time_t a_date = 3600 * 24 * (31 + 4); // Feb 5th
    time_t a_time = 3600 * 15 + 60 * 33;  // 15:33:13
    time_t a_timesec = 13;
    time_t a_datetime = a_date + a_time + a_timesec;

    std::string pat[] = {"a",        "Thu", "A",  "Thursday", "b",   "Feb", "B",    "February",    "d",  "05",    "D",
                         "02/05/70", "e",   "5",  "h",        "Feb", "H",   "15",   "I",           "03", "m",     "02",
                         "M",        "33",  "n",  "\n",       "p",   "PM",  "r",    "03:33:13 PM", "R",  "15:33", "S",
                         "13",       "t",   "\t", "y",        "70",  "Y",   "1970", "%",           "%"};

    for(unsigned i = 0; i < sizeof(pat) / sizeof(pat[0]); i += 2) {
        ss.str("");
        ss << boost::locale::as::ftime("%" + pat[i]) << a_datetime;
        TEST_EQ(ss.str(), pat[i + 1]);
    }
}

void test_main(int /*argc*/, char** /*argv*/)
{
#ifdef BOOST_LOCALE_NO_WINAPI_BACKEND
    std::cout << "WinAPI Backend is not build... Skipping\n";
    return;
#endif
    boost::locale::localization_backend_manager mgr = boost::locale::localization_backend_manager::global();
    mgr.select("winapi");
    boost::locale::localization_backend_manager::global(mgr);
    boost::locale::generator gen;

    for(const auto& name_lcid : {std::make_pair("en_US.UTF-8", 0x0409),
                                 std::make_pair("he_IL.UTF-8", 0x040D),
                                 std::make_pair("ru_RU.UTF-8", 0x0419)})
    {
        const std::string name = name_lcid.first;
        std::cout << "- " << name << " locale" << std::endl;
        if(boost::locale::impl_win::locale_to_lcid(name) == 0) {
            std::cout << "-- not supported, skipping" << std::endl;
            continue;
        }
        std::locale l1 = gen(name);
        std::cout << "-- UTF-8" << std::endl;
        test_by_char<char>(l1, name, name_lcid.second);
        std::cout << "-- UTF-16" << std::endl;
        test_by_char<wchar_t>(l1, name, name_lcid.second);
    }
    std::cout << "- Testing strftime" << std::endl;
    test_date_time(gen("en_US.UTF-8"));
}

// boostinspect:noascii
