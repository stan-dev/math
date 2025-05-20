//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale/format.hpp>
#include <boost/locale/formatting.hpp>
#include <boost/locale/generator.hpp>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "boostLocale/test/tools.hpp"
#include "boostLocale/test/unit_test.hpp"

BOOST_LOCALE_DISABLE_UNREACHABLE_CODE_WARNING
void test_main(int /*argc*/, char** /*argv*/)
{
#ifndef BOOST_LOCALE_WITH_ICU
    std::cout << "ICU is not build... Skipping\n";
    return;
#endif
    const time_t now = std::time(nullptr);
    boost::locale::generator gen;
    std::locale::global(gen("en_US.UTF-8"));

    for(int i = 0; i < 366; i++) {
        time_t point = now + i * 24 * 3600;
        std::stringstream ss;
        ss << boost::locale::format("{1,ftime='%H %M %S'}") % point;
        int icu_hour = 0, icu_min = 0, icu_sec = 0;
        ss >> icu_hour >> icu_min >> icu_sec;
        std::tm* tm = localtime_wrap(&point);
        TEST_EQ(icu_hour, tm->tm_hour);
        TEST_EQ(icu_min, tm->tm_min);
        TEST_EQ(icu_sec, tm->tm_sec);
    }
}

// boostinspect:noascii
