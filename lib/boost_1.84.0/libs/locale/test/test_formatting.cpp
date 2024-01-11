//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
// Copyright (c) 2021-2022 Alexander Grund
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale/date_time.hpp>
#include <boost/locale/encoding_utf.hpp>
#include <boost/locale/format.hpp>
#include <boost/locale/formatting.hpp>
#include <boost/locale/generator.hpp>
#include <boost/locale/localization_backend.hpp>
#include <cstdint>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>

#include "boostLocale/test/tools.hpp"
#include "boostLocale/test/unit_test.hpp"

const std::string test_locale_name = "en_US";
std::string message_path = "./";

#ifdef BOOST_LOCALE_WITH_ICU
#    include <unicode/datefmt.h>
#    include <unicode/numfmt.h>
#    include <unicode/timezone.h>
#    include <unicode/uversion.h>
#    define BOOST_LOCALE_ICU_VERSION (U_ICU_VERSION_MAJOR_NUM * 100 + U_ICU_VERSION_MINOR_NUM)
#    define BOOST_LOCALE_ICU_VERSION_EXACT (BOOST_LOCALE_ICU_VERSION * 100 + U_ICU_VERSION_PATCHLEVEL_NUM)

const icu::Locale& get_icu_test_locale()
{
    static icu::Locale locale = icu::Locale::createCanonical(test_locale_name.c_str());
    return locale;
}

std::string from_icu_string(const icu::UnicodeString& str)
{
    return boost::locale::conv::utf_to_utf<char>(str.getBuffer(), str.getBuffer() + str.length());
}
#else
#    define BOOST_LOCALE_ICU_VERSION 0
#    define BOOST_LOCALE_ICU_VERSION_EXACT 0
#endif

// Currency style changes between ICU versions, so get "real" value from ICU
#if BOOST_LOCALE_ICU_VERSION >= 402

std::string get_icu_currency_iso(const double value)
{
#    if BOOST_LOCALE_ICU_VERSION >= 408
    auto styleIso = UNUM_CURRENCY_ISO;
#    else
    auto styleIso = icu::NumberFormat::kIsoCurrencyStyle;
#    endif
    UErrorCode err = U_ZERO_ERROR;
    std::unique_ptr<icu::NumberFormat> fmt(icu::NumberFormat::createInstance(get_icu_test_locale(), styleIso, err));
    TEST_REQUIRE(U_SUCCESS(err) && fmt.get());

    icu::UnicodeString tmp;
    return from_icu_string(fmt->format(value, tmp));
}

#endif

using format_style_t = std::ios_base&(std::ios_base&);

#ifdef BOOST_LOCALE_WITH_ICU
std::string get_icu_gmt_name(icu::TimeZone::EDisplayType style)
{
    icu::UnicodeString tmp;
    return from_icu_string(icu::TimeZone::getGMT()->getDisplayName(false, style, get_icu_test_locale(), tmp));
}

// This changes between ICU versions, e.g. "GMT" or "Greenwich Mean Time"
const std::string icu_full_gmt_name = get_icu_gmt_name(icu::TimeZone::EDisplayType::LONG);

std::string get_ICU_time(format_style_t style, const time_t ts, const char* tz = nullptr)
{
    using icu::DateFormat;
    DateFormat::EStyle icu_style = DateFormat::kDefault;
    namespace as = boost::locale::as;
    if(style == as::time_short)
        icu_style = DateFormat::kShort;
    else if(style == as::time_medium)
        icu_style = DateFormat::kMedium;
    else if(style == as::time_long)
        icu_style = DateFormat::kLong;
    else if(style == as::time_full)
        icu_style = DateFormat::kFull;
    std::unique_ptr<icu::DateFormat> fmt(icu::DateFormat::createTimeInstance(icu_style, get_icu_test_locale()));
    if(!tz)
        fmt->setTimeZone(*icu::TimeZone::getGMT());
    else
        fmt->adoptTimeZone(icu::TimeZone::createTimeZone(icu::UnicodeString::fromUTF8(tz)));
    icu::UnicodeString s;
    return from_icu_string(fmt->format(ts * 1000., s));
}

std::string get_ICU_date(format_style_t style, const time_t ts)
{
    using icu::DateFormat;
    DateFormat::EStyle icu_style = DateFormat::kDefault;
    namespace as = boost::locale::as;
    if(style == as::date_short)
        icu_style = DateFormat::kShort;
    else if(style == as::date_medium)
        icu_style = DateFormat::kMedium;
    else if(style == as::date_long)
        icu_style = DateFormat::kLong;
    else if(style == as::date_full)
        icu_style = DateFormat::kFull;
    std::unique_ptr<icu::DateFormat> fmt(icu::DateFormat::createDateInstance(icu_style, get_icu_test_locale()));
    fmt->setTimeZone(*icu::TimeZone::getGMT());
    icu::UnicodeString s;
    return from_icu_string(fmt->format(ts * 1000., s));
}

std::string get_ICU_datetime(format_style_t style, const time_t ts)
{
    using icu::DateFormat;
    DateFormat::EStyle icu_style = DateFormat::kDefault;
    namespace as = boost::locale::as;
    if(style == as::time_short)
        icu_style = DateFormat::kShort;
    else if(style == as::time_medium)
        icu_style = DateFormat::kMedium;
    else if(style == as::time_long)
        icu_style = DateFormat::kLong;
    else if(style == as::time_full)
        icu_style = DateFormat::kFull;
    std::unique_ptr<icu::DateFormat> fmt(
      icu::DateFormat::createDateTimeInstance(icu_style, icu_style, get_icu_test_locale()));
    fmt->setTimeZone(*icu::TimeZone::getGMT());
    icu::UnicodeString s;
    return from_icu_string(fmt->format(ts * 1000., s));
}

#else
const std::string icu_full_gmt_name;
// clang-format off
std::string get_ICU_time(...){ return ""; } // LCOV_EXCL_LINE
std::string get_ICU_datetime(...){ return ""; } // LCOV_EXCL_LINE
std::string get_ICU_date(...){ return ""; } // LCOV_EXCL_LINE
// clang-format on
#endif

using namespace boost::locale;

template<typename CharType, typename T>
void test_fmt_impl(std::basic_ostringstream<CharType>& ss,
                   const T& value,
                   const std::basic_string<CharType>& expected,
                   int line)
{
    ss << value;
    test_eq_impl(ss.str(), expected, "", line);
}

template<typename T, typename CharType>
void test_parse_impl(std::basic_istringstream<CharType>& ss, const T& expected, int line)
{
    T v;
    ss >> v >> std::ws;
    test_eq_impl(v, expected, "v == expected", line);
    test_eq_impl(ss.eof(), true, "ss.eof()", line);
}

template<typename T, typename CharType>
void test_parse_at_impl(std::basic_istringstream<CharType>& ss, const T& expected, int line)
{
    T v;
    CharType c_at;
    ss >> v >> std::skipws >> c_at;
    test_eq_impl(v, expected, "v == expected", line);
    test_eq_impl(c_at, '@', "c_at == @", line);
}

template<typename T, typename CharType>
void test_parse_fail_impl(std::basic_istringstream<CharType>& ss, int line)
{
    T v;
    ss >> v;
    test_eq_impl(ss.fail(), true, "ss.fail()", line);
}

#define TEST_FMT(manip, value, expected)                                                  \
    do {                                                                                  \
        std::basic_ostringstream<CharType> ss;                                            \
        ss.imbue(loc);                                                                    \
        ss << manip;                                                                      \
        test_fmt_impl(ss, (value), to_correct_string<CharType>(expected, loc), __LINE__); \
        BOOST_LOCALE_START_CONST_CONDITION                                                \
    } while(0) BOOST_LOCALE_END_CONST_CONDITION

#define TEST_PARSE_FAILS(manip, actual, type)             \
    do {                                                  \
        std::basic_istringstream<CharType> ss;            \
        ss.imbue(loc);                                    \
        ss.str(to_correct_string<CharType>(actual, loc)); \
        ss >> manip;                                      \
        test_parse_fail_impl<type>(ss, __LINE__);         \
        BOOST_LOCALE_START_CONST_CONDITION                \
    } while(0) BOOST_LOCALE_END_CONST_CONDITION

#define TEST_PARSE(manip, value, expected)                              \
    do {                                                                \
        const auto str_value = to_correct_string<CharType>(value, loc); \
        {                                                               \
            std::basic_istringstream<CharType> ss;                      \
            ss.imbue(loc);                                              \
            ss.str(str_value);                                          \
            ss >> manip;                                                \
            test_parse_impl(ss, expected, __LINE__);                    \
        }                                                               \
        {                                                               \
            std::basic_istringstream<CharType> ss;                      \
            ss.imbue(loc);                                              \
            ss.str(str_value + CharType('@'));                          \
            ss >> manip;                                                \
            test_parse_at_impl(ss, expected, __LINE__);                 \
        }                                                               \
        BOOST_LOCALE_START_CONST_CONDITION                              \
    } while(0) BOOST_LOCALE_END_CONST_CONDITION

#define TEST_FMT_PARSE_1(manip, value_in, value_str) \
    do {                                             \
        const std::string value_str_ = value_str;    \
        TEST_FMT(manip, value_in, value_str_);       \
        TEST_PARSE(manip, value_str_, value_in);     \
        BOOST_LOCALE_START_CONST_CONDITION           \
    } while(0) BOOST_LOCALE_END_CONST_CONDITION

#define TEST_FMT_PARSE_2(m1, m2, value_in, value_str) \
    do {                                              \
        const std::string value_str_ = value_str;     \
        TEST_FMT(m1 << m2, value_in, value_str_);     \
        TEST_PARSE(m1 >> m2, value_str_, value_in);   \
        BOOST_LOCALE_START_CONST_CONDITION            \
    } while(0) BOOST_LOCALE_END_CONST_CONDITION

#define TEST_FMT_PARSE_2_2(m1, m2, value_in, value_str, value_parsed) \
    do {                                                              \
        const std::string value_str_ = value_str;                     \
        TEST_FMT(m1 << m2, value_in, value_str_);                     \
        TEST_PARSE(m1 >> m2, value_str_, value_parsed);               \
        BOOST_LOCALE_START_CONST_CONDITION                            \
    } while(0) BOOST_LOCALE_END_CONST_CONDITION

#define TEST_FMT_PARSE_3(m1, m2, m3, value_in, value_str) \
    do {                                                  \
        const std::string value_str_ = value_str;         \
        TEST_FMT(m1 << m2 << m3, value_in, value_str_);   \
        TEST_PARSE(m1 >> m2 >> m3, value_str_, value_in); \
        BOOST_LOCALE_START_CONST_CONDITION                \
    } while(0) BOOST_LOCALE_END_CONST_CONDITION

#define TEST_FMT_PARSE_3_2(m1, m2, m3, value_in, value_str, value_parsed) \
    do {                                                                  \
        const std::string value_str_ = value_str;                         \
        TEST_FMT(m1 << m2 << m3, value_in, value_str_);                   \
        TEST_PARSE(m1 >> m2 >> m3, value_str_, value_parsed);             \
        BOOST_LOCALE_START_CONST_CONDITION                                \
    } while(0) BOOST_LOCALE_END_CONST_CONDITION

#define TEST_FMT_PARSE_4(m1, m2, m3, m4, value_in, value_str)   \
    do {                                                        \
        const std::string value_str_ = value_str;               \
        TEST_FMT(m1 << m2 << m3 << m4, value_in, value_str_);   \
        TEST_PARSE(m1 >> m2 >> m3 >> m4, value_str_, value_in); \
        BOOST_LOCALE_START_CONST_CONDITION                      \
    } while(0) BOOST_LOCALE_END_CONST_CONDITION

#define TEST_FMT_PARSE_4_2(m1, m2, m3, m4, value_in, value_str, value_parsed) \
    do {                                                                      \
        const std::string value_str_ = value_str;                             \
        TEST_FMT(m1 << m2 << m3 << m4, value_in, value_str_);                 \
        TEST_PARSE(m1 >> m2 >> m3 >> m4, value_str_, value_parsed);           \
        BOOST_LOCALE_START_CONST_CONDITION                                    \
    } while(0) BOOST_LOCALE_END_CONST_CONDITION

#define TEST_MIN_MAX_FMT(as, type, minval, maxval)          \
    TEST_FMT(as, std::numeric_limits<type>::min(), minval); \
    TEST_FMT(as, std::numeric_limits<type>::max(), maxval)

#define TEST_MIN_MAX_PARSE(as, type, minval, maxval)          \
    TEST_PARSE(as, minval, std::numeric_limits<type>::min()); \
    TEST_PARSE(as, maxval, std::numeric_limits<type>::max())

#define TEST_MIN_MAX(type, minval, maxval)              \
    TEST_MIN_MAX_FMT(as::number, type, minval, maxval); \
    TEST_MIN_MAX_PARSE(as::number, type, minval, maxval)

#define TEST_MIN_MAX_POSIX(type)                                                      \
    do {                                                                              \
        const std::string minval = as_posix_string(std::numeric_limits<type>::min()); \
        const std::string maxval = as_posix_string(std::numeric_limits<type>::max()); \
        TEST_MIN_MAX_FMT(as::posix, type, minval, maxval);                            \
        TEST_MIN_MAX_PARSE(as::posix, type, minval, maxval);                          \
        BOOST_LOCALE_START_CONST_CONDITION                                            \
    } while(0) BOOST_LOCALE_END_CONST_CONDITION

bool stdlib_correctly_errors_on_out_of_range_int16()
{
    static bool fails = []() -> bool {
        std::stringstream ss("65000");
        ss.imbue(std::locale::classic());
        int16_t v = 0;
        ss >> v;
        return ss.fail();
    }();
    return fails;
}

template<typename T>
std::string as_posix_string(const T v)
{
    std::ostringstream ss;
    ss.imbue(std::locale::classic());
    ss << v;
    return ss.str();
}

template<typename CharType>
void test_as_posix(const std::string& e_charset = "UTF-8")
{
    using boost::locale::localization_backend_manager;
    const auto orig_backend = localization_backend_manager::global();
    for(const std::string& backendName : orig_backend.get_all_backends()) {
        std::cout << "Backend: " << backendName << std::endl;
        auto backend = orig_backend;
        backend.select(backendName);
        localization_backend_manager::global(backend);
        for(const std::string name : {"en_US", "ru_RU", "de_DE"}) {
            const std::locale loc = boost::locale::generator{}(name + "." + e_charset);
            TEST_MIN_MAX_POSIX(int16_t);
            TEST_MIN_MAX_POSIX(uint16_t);

            TEST_MIN_MAX_POSIX(int32_t);
            TEST_MIN_MAX_POSIX(uint32_t);
            TEST_MIN_MAX_POSIX(signed long);
            TEST_MIN_MAX_POSIX(unsigned long);

            TEST_MIN_MAX_POSIX(int64_t);
            TEST_MIN_MAX_POSIX(uint64_t);
            TEST_MIN_MAX_POSIX(signed long long);
            TEST_MIN_MAX_POSIX(unsigned long long);

            TEST_FMT_PARSE_1(as::posix, 1.25f, "1.25");
            TEST_FMT_PARSE_1(as::posix, -4.57, "-4.57");
            TEST_FMT_PARSE_1(as::posix, 3.815l, "3.815");
        }
    }
    localization_backend_manager::global(orig_backend);
}

template<typename CharType>
void test_manip(std::string e_charset = "UTF-8")
{
    test_as_posix<CharType>(e_charset);
    using string_type = std::basic_string<CharType>;
    boost::locale::generator g;
    for(const auto& name_number : {std::make_pair("en_US", "1,200.1"),
                                   std::make_pair("he_IL", "1,200.1"),
                                   std::make_pair("ru_RU",
                                                  "1\xC2\xA0"
                                                  "200,1")})
    {
        const std::string locName = std::string(name_number.first) + "." + e_charset;
        std::cout << "-- " << locName << '\n';
        const std::locale loc = g(locName);
        TEST_FMT_PARSE_1(as::posix, 1200.1, "1200.1");
        TEST_FMT_PARSE_1(as::number, 1200.1, name_number.second);
    }

    const std::locale loc = g(test_locale_name + "." + e_charset);
    TEST_FMT_PARSE_1(as::posix, 1200.1, "1200.1");
    TEST_FMT_PARSE_1(as::number, 1200.1, "1,200.1");
    TEST_FMT(as::number << std::setfill(CharType('_')) << std::setw(6), 1534, "_1,534");
    TEST_FMT(as::number << std::left << std::setfill(CharType('_')) << std::setw(6), 1534, "1,534_");

    // Ranges
    TEST_MIN_MAX(int16_t, "-32,768", "32,767");
    TEST_MIN_MAX(uint16_t, "0", "65,535");
    TEST_PARSE_FAILS(as::number, "-1", uint16_t);
    if(stdlib_correctly_errors_on_out_of_range_int16())
        TEST_PARSE_FAILS(as::number, "65,535", int16_t);

    TEST_MIN_MAX(int32_t, "-2,147,483,648", "2,147,483,647");
    TEST_MIN_MAX(uint32_t, "0", "4,294,967,295");
    TEST_PARSE_FAILS(as::number, "-1", uint32_t);
    TEST_PARSE_FAILS(as::number, "4,294,967,295", int32_t);

    TEST_MIN_MAX(int64_t, "-9,223,372,036,854,775,808", "9,223,372,036,854,775,807");
    // ICU does not support uint64, but we have a fallback to format it at least
    TEST_MIN_MAX_FMT(as::number, uint64_t, "0", "18446744073709551615");
    TEST_PARSE_FAILS(as::number, "-1", uint64_t);

    TEST_FMT_PARSE_3(as::number, std::left, std::setw(3), 15, "15 ");
    TEST_FMT_PARSE_3(as::number, std::right, std::setw(3), 15, " 15");
    TEST_FMT_PARSE_3(as::number, std::setprecision(3), std::fixed, 13.1, "13.100");
    TEST_FMT_PARSE_3(as::number, std::setprecision(3), std::scientific, 13.1, "1.310E1");

    TEST_PARSE_FAILS(as::number, "", int);
    TEST_PARSE_FAILS(as::number, "--3", int);
    TEST_PARSE_FAILS(as::number, "y", int);

    TEST_FMT_PARSE_1(as::percent, 0.1, "10%");
    TEST_FMT_PARSE_3(as::percent, std::fixed, std::setprecision(1), 0.10, "10.0%");

    TEST_PARSE_FAILS(as::percent, "1", double);

    TEST_FMT_PARSE_1(as::currency, 1345, "$1,345.00");
    TEST_FMT_PARSE_1(as::currency, 1345.34, "$1,345.34");

    TEST_PARSE_FAILS(as::currency, "$", double);

#if BOOST_LOCALE_ICU_VERSION >= 402
    TEST_FMT_PARSE_2(as::currency, as::currency_national, 1345, "$1,345.00");
    TEST_FMT_PARSE_2(as::currency, as::currency_national, 1345.34, "$1,345.34");
    TEST_FMT_PARSE_2(as::currency, as::currency_iso, 1345, get_icu_currency_iso(1345));
    TEST_FMT_PARSE_2(as::currency, as::currency_iso, 1345.34, get_icu_currency_iso(1345.34));
#endif
    TEST_FMT_PARSE_1(as::spellout, 10, "ten");
#if 402 <= BOOST_LOCALE_ICU_VERSION && BOOST_LOCALE_ICU_VERSION < 408
    if(e_charset == "UTF-8")
        TEST_FMT(as::ordinal, 1, "1\xcb\xa2\xe1\xb5\x97"); // 1st with st as ligatures
#else
    TEST_FMT(as::ordinal, 1, "1st");
#endif

    time_t a_date = 3600 * 24 * (31 + 4); // Feb 5th
    time_t a_time = 3600 * 15 + 60 * 33;  // 15:33:05
    time_t a_timesec = 13;
    time_t a_datetime = a_date + a_time + a_timesec;

    TEST_FMT_PARSE_2_2(as::date, as::gmt, a_datetime, "Feb 5, 1970", a_date);
    TEST_FMT_PARSE_3_2(as::date, as::date_short, as::gmt, a_datetime, "2/5/70", a_date);
    TEST_FMT_PARSE_3_2(as::date, as::date_medium, as::gmt, a_datetime, "Feb 5, 1970", a_date);
    TEST_FMT_PARSE_3_2(as::date, as::date_long, as::gmt, a_datetime, "February 5, 1970", a_date);
    TEST_FMT_PARSE_3_2(as::date, as::date_full, as::gmt, a_datetime, "Thursday, February 5, 1970", a_date);

    TEST_PARSE_FAILS(as::date >> as::date_short, "aa/bb/cc", double);

    std::string icu_time_def = get_ICU_time(as::time, a_datetime);
    std::string icu_time_short = get_ICU_time(as::time_short, a_datetime);
    std::string icu_time_medium = get_ICU_time(as::time_medium, a_datetime);
    std::string icu_time_long = get_ICU_time(as::time_long, a_datetime);
    std::string icu_time_full = get_ICU_time(as::time_full, a_datetime);

    TEST_PARSE(as::time >> as::gmt, "3:33:13 PM", a_time + a_timesec);
    TEST_FMT_PARSE_2_2(as::time, as::gmt, a_datetime, icu_time_def, a_time + a_timesec);

    TEST_PARSE(as::time >> as::time_short >> as::gmt, "3:33 PM", a_time);
    TEST_FMT_PARSE_3_2(as::time, as::time_short, as::gmt, a_datetime, icu_time_short, a_time);

    TEST_PARSE(as::time >> as::time_medium >> as::gmt, "3:33:13 PM", a_time + a_timesec);
    TEST_FMT_PARSE_3_2(as::time, as::time_medium, as::gmt, a_datetime, icu_time_medium, a_time + a_timesec);

    TEST_PARSE(as::time >> as::time_long >> as::gmt, "3:33:13 PM GMT", a_time + a_timesec);
    TEST_FMT_PARSE_3_2(as::time, as::time_long, as::gmt, a_datetime, icu_time_long, a_time + a_timesec);
    // ICU 4.8.0 has a bug which makes parsing the full time fail when anything follows the time zone
#if BOOST_LOCALE_ICU_VERSION_EXACT != 40800
    TEST_PARSE(as::time >> as::time_full >> as::gmt, "3:33:13 PM GMT+00:00", a_time + a_timesec);
    TEST_FMT_PARSE_3_2(as::time, as::time_full, as::gmt, a_datetime, icu_time_full, a_time + a_timesec);
#endif
    TEST_PARSE_FAILS(as::time, "AM", double);

    icu_time_def = get_ICU_time(as::time, a_datetime, "GMT+01:00");
    icu_time_short = get_ICU_time(as::time_short, a_datetime, "GMT+01:00");
    icu_time_medium = get_ICU_time(as::time_medium, a_datetime, "GMT+01:00");
    icu_time_long = get_ICU_time(as::time_long, a_datetime, "GMT+01:00");
    icu_time_full = get_ICU_time(as::time_full, a_datetime, "GMT+01:00");

    TEST_PARSE(as::time >> as::time_zone("GMT+01:00"), "4:33:13 PM", a_time + a_timesec);
    TEST_FMT_PARSE_2_2(as::time, as::time_zone("GMT+01:00"), a_datetime, icu_time_def, a_time + a_timesec);

    TEST_PARSE(as::time >> as::time_short >> as::time_zone("GMT+01:00"), "4:33 PM", a_time);
    TEST_FMT_PARSE_3_2(as::time, as::time_short, as::time_zone("GMT+01:00"), a_datetime, icu_time_short, a_time);

    TEST_PARSE(as::time >> as::time_medium >> as::time_zone("GMT+01:00"), "4:33:13 PM", a_time + a_timesec);
    TEST_FMT_PARSE_3_2(as::time,
                       as::time_medium,
                       as::time_zone("GMT+01:00"),
                       a_datetime,
                       icu_time_medium,
                       a_time + a_timesec);

    TEST_PARSE(as::time >> as::time_long >> as::time_zone("GMT+01:00"), "4:33:13 PM GMT+01:00", a_time + a_timesec);
    TEST_FMT_PARSE_3_2(as::time,
                       as::time_long,
                       as::time_zone("GMT+01:00"),
                       a_datetime,
                       icu_time_long,
                       a_time + a_timesec);
#if !(BOOST_LOCALE_ICU_VERSION == 308 && defined(__CYGWIN__)) // Known failure due to ICU issue
    TEST_PARSE(as::time >> as::time_full >> as::time_zone("GMT+01:00"), "4:33:13 PM GMT+01:00", a_time + a_timesec);
    TEST_FMT_PARSE_3_2(as::time,
                       as::time_full,
                       as::time_zone("GMT+01:00"),
                       a_datetime,
                       icu_time_full,
                       a_time + a_timesec);
#endif

    const std::string icu_def = get_ICU_datetime(as::time, a_datetime);
    const std::string icu_short = get_ICU_datetime(as::time_short, a_datetime);
    const std::string icu_medium = get_ICU_datetime(as::time_medium, a_datetime);
    const std::string icu_long = get_ICU_datetime(as::time_long, a_datetime);
    const std::string icu_full = get_ICU_datetime(as::time_full, a_datetime);

    TEST_PARSE(as::datetime >> as::gmt, "Feb 5, 1970 3:33:13 PM", a_datetime);
    TEST_FMT_PARSE_2(as::datetime, as::gmt, a_datetime, icu_def);

    TEST_PARSE(as::datetime >> as::date_short >> as::time_short >> as::gmt, "2/5/70 3:33 PM", a_date + a_time);
    TEST_FMT_PARSE_4_2(as::datetime, as::date_short, as::time_short, as::gmt, a_datetime, icu_short, a_date + a_time);

    TEST_PARSE(as::datetime >> as::date_medium >> as::time_medium >> as::gmt, "Feb 5, 1970 3:33:13 PM", a_datetime);
    TEST_FMT_PARSE_4(as::datetime, as::date_medium, as::time_medium, as::gmt, a_datetime, icu_medium);

    TEST_PARSE(as::datetime >> as::date_long >> as::time_long >> as::gmt,
               "February 5, 1970 3:33:13 PM GMT",
               a_datetime);
    TEST_FMT_PARSE_4(as::datetime, as::date_long, as::time_long, as::gmt, a_datetime, icu_long);
#if BOOST_LOCALE_ICU_VERSION_EXACT != 40800
    // ICU 4.8.0 has a bug which makes parsing the full time fail when anything follows the time zone
    TEST_PARSE(as::datetime >> as::date_full >> as::time_full >> as::gmt,
               "Thursday, February 5, 1970 3:33:13 PM Greenwich Mean Time",
               a_datetime);
    TEST_FMT_PARSE_4(as::datetime, as::date_full, as::time_full, as::gmt, a_datetime, icu_full);
#endif

    const std::pair<char, std::string> mark_test_cases[] = {
      std::make_pair('a', "Thu"),
      std::make_pair('A', "Thursday"),
      std::make_pair('b', "Feb"),
      std::make_pair('B', "February"),
      std::make_pair('c', icu_full),
      std::make_pair('d', "05"),
      std::make_pair('e', "5"),
      std::make_pair('h', "Feb"),
      std::make_pair('H', "15"),
      std::make_pair('I', "03"),
      std::make_pair('j', "36"),
      std::make_pair('m', "02"),
      std::make_pair('M', "33"),
      std::make_pair('n', "\n"),
      std::make_pair('p', "PM"),
      std::make_pair('r', "03:33:13 PM"),
      std::make_pair('R', "15:33"),
      std::make_pair('S', "13"),
      std::make_pair('t', "\t"),
      std::make_pair('T', "15:33:13"),
      std::make_pair('x', "Feb 5, 1970"),
      std::make_pair('X', get_ICU_time(as::time, a_datetime)),
      std::make_pair('y', "70"),
      std::make_pair('Y', "1970"),
      std::make_pair('Z', icu_full_gmt_name),
      std::make_pair('%', "%"),
    };

    for(const auto& mark_result : mark_test_cases) {
        string_type format_string;
        format_string += static_cast<CharType>('%');
        format_string += static_cast<CharType>(mark_result.first);
        std::cout << "Test: %" << mark_result.first << "\n";
        std::basic_ostringstream<CharType> ss;
        ss.imbue(loc);
        ss << as::ftime(format_string) << as::gmt << a_datetime;
        TEST_EQ(ss.str(), to_correct_string<CharType>(mark_result.second, loc));
    }

    {
        const time_t now = time(nullptr);
        boost::locale::time_zone::global("GMT+4:00");
        const time_t local_now = now + 3600 * 4;
        char time_str[256];
        const std::string format = "%H:%M:%S";
        const string_type format_string(format.begin(), format.end());

        std::basic_ostringstream<CharType> ss;
        ss.imbue(loc);
        // By default the globally set (local) time zone is used
        ss << as::ftime(format_string) << now;
        strftime(time_str, sizeof(time_str), format.c_str(), gmtime_wrap(&local_now));
        TEST_EQ(ss.str(), to<CharType>(time_str));
        // We can manually tell it to use the local time zone
        empty_stream(ss) << as::ftime(format_string) << as::local_time << now;
        TEST_EQ(ss.str(), to<CharType>(time_str));
        // Or e.g. GMT
        empty_stream(ss) << as::ftime(format_string) << as::gmt << now;
        strftime(time_str, sizeof(time_str), format.c_str(), gmtime_wrap(&now));
        TEST_EQ(ss.str(), to<CharType>(time_str));
    }
    const std::pair<std::string, std::string> format_string_test_cases[] = {
      std::make_pair("Now is %A, %H o'clo''ck ' or not ' ", "Now is Thursday, 15 o'clo''ck ' or not ' "),
      std::make_pair("'test %H'", "'test 15'"),
      std::make_pair("%H'", "15'"),
      std::make_pair("'%H'", "'15'"),
    };

    for(const auto& test_case : format_string_test_cases) {
        const string_type format_string(test_case.first.begin(), test_case.first.end());
        std::cout << "Test: '" << test_case.first << "'\n";
        std::basic_ostringstream<CharType> ss;
        ss.imbue(loc);
        ss << as::ftime(format_string) << as::gmt << a_datetime;
        TEST_EQ(ss.str(), to<CharType>(test_case.second));
    }
}

template<typename CharType, typename... Ts>
std::basic_string<CharType> do_format(const std::locale& loc, const std::basic_string<CharType> fmt_str, Ts&&... ts)
{
    boost::locale::basic_format<CharType> fmt(fmt_str);
    using expander = int[];
    (void)expander{0, (fmt % std::forward<Ts>(ts), 0)...};
    return fmt.str(loc);
}

template<typename CharType, size_t size, typename... Ts>
std::basic_string<CharType> do_format(const std::locale& loc, const char (&fmt_str)[size], Ts&&... ts)
{
    return do_format(loc, ascii_to<CharType>(fmt_str), std::forward<Ts>(ts)...);
}

template<typename CharType, typename T>
void test_format_class_impl(const std::string& fmt_string,
                            const T& value,
                            const std::string& expected_str,
                            const std::locale& loc,
                            unsigned line)
{
    using format_type = boost::locale::basic_format<CharType>;
    format_type fmt(std::basic_string<CharType>(fmt_string.begin(), fmt_string.end()));
    fmt % value;
    std::basic_string<CharType> expected_str_loc(to_correct_string<CharType>(expected_str, loc));
    test_eq_impl(fmt.str(loc), expected_str_loc, ("Format: " + fmt_string).c_str(), line);
}

template<typename CharType>
void test_format_class(std::string charset = "UTF-8")
{
    using string_type = std::basic_string<CharType>;
    using format_type = boost::locale::basic_format<CharType>;

    boost::locale::generator g;
    std::locale loc = g(test_locale_name + "." + charset);

    // Simple tests using same input/output
    {
        const string_type fmt_string = ascii_to<CharType>("{3} {1} {2}");
        const string_type expected = ascii_to<CharType>("3 1 2");

        // Output format to stream
        {
            std::basic_ostringstream<CharType> ss;
            ss.imbue(loc);

            // Stream formatted output
            empty_stream(ss) << format_type(fmt_string) % 1 % 2 % 3;
            TEST_EQ(ss.str(), expected);

            // Stream translated output
            empty_stream(ss) << format_type(boost::locale::translate(fmt_string)) % 1 % 2 % 3;
            TEST_EQ(ss.str(), expected);
        }

        // Multi-step: Create, format & output via str() method
        {
            format_type fmt(fmt_string);
            TEST_EQ((fmt % 1 % 2 % 3).str(loc), expected);
        }
        // Output via str() on intermediate with ctor from string and C-string
        TEST_EQ((format_type(fmt_string.c_str()) % 1 % 2 % 3).str(loc), expected);
        TEST_EQ((format_type(fmt_string) % 1 % 2 % 3).str(loc), expected);
    }

    // Actually translate something
    {
        g.add_messages_domain("default/ISO-8859-8");
        g.add_messages_path(message_path);
        std::locale loc_he = g("he_IL.UTF-8");
        const string_type hello = ascii_to<CharType>("hello");
        const string_type hello_he = to_correct_string<CharType>("שלום", loc_he);
        format_type fmt(boost::locale::translate(hello));
        TEST_EQ(fmt.str(loc_he), hello_he);
        // Use current global locale if none is given to str()
        std::locale old_locale = std::locale::global(g("en_US.UTF-8"));
        TEST_EQ(fmt.str(), hello); // Not translated in en_US
        std::locale::global(loc_he);
        TEST_EQ(fmt.str(), hello_he); // translated in he_IL

        // Movable
        {
            format_type fmt2 = format_type(ascii_to<CharType>("{3} {1} {2}"));
            const int i1 = 1, i2 = 2, i3 = 3, i42 = 42;
            fmt2 % i1 % i2 % i3;
            fmt2 = format_type(ascii_to<CharType>("{1}"));
            TEST_EQ(fmt2.str(), ascii_to<CharType>("")); // No bound value
            TEST_EQ((fmt2 % i42).str(), ascii_to<CharType>("42"));
            // Can't move with bound params
            TEST_THROWS(format_type fmt3(std::move(fmt2)), std::exception);
            TEST_EQ(fmt2.str(), ascii_to<CharType>("42")); // Original unchanged
            fmt2 = format_type(ascii_to<CharType>("{1}"));
            fmt2 % i1;
            format_type fmt3(string_type{});
            TEST_THROWS(fmt3 = std::move(fmt2), std::exception);
            fmt2 = format_type(ascii_to<CharType>("{1}"));
            fmt3 = std::move(fmt2);
            TEST_EQ((fmt3 % 42).str(), ascii_to<CharType>("42"));

            fmt2 = format_type(hello);
            TEST_EQ(fmt2.str(), hello); // Not translated
            fmt2 = format_type(boost::locale::translate(hello));
            TEST_EQ(fmt2.str(), hello_he); // Translated
        }
        // Restore
        std::locale::global(old_locale);
    }
    // Allows many params
    TEST_EQ(do_format<CharType>(loc, "{1}{2}{3}{4}{5}{10}{9}{8}{7}{6}", 11, 22, 33, 44, 55, 'a', 'b', 'c', 'd', 'f'),
            ascii_to<CharType>("1122334455fdcba"));
    // Not passed placeholders are removed
    TEST_EQ(do_format<CharType>(loc, "{1}{3}{2}", "hello", "world"), ascii_to<CharType>("helloworld"));
    TEST_EQ(do_format<CharType>(loc, "{1}"), ascii_to<CharType>(""));
    // Invalid indices are ignored
    TEST_EQ(do_format<CharType>(loc, "b{}e"), ascii_to<CharType>("be"));
    TEST_EQ(do_format<CharType>(loc, "b{0}e", 1), ascii_to<CharType>("be"));
    TEST_EQ(do_format<CharType>(loc, "b{-1}e", 1), ascii_to<CharType>("be"));
    TEST_EQ(do_format<CharType>(loc, "b{1.x}e"), ascii_to<CharType>("be"));
    // Unexpected closing brace and other chars are ignored
    TEST_EQ(do_format<CharType>(loc, " = , } 3"), ascii_to<CharType>(" = , } 3"));
    // Trailing opening brace is ignored
    TEST_EQ(do_format<CharType>(loc, "End {"), ascii_to<CharType>("End "));
    // Trailing closing brace is added like any other char
    TEST_EQ(do_format<CharType>(loc, "End}"), ascii_to<CharType>("End}"));
    // Escaped trailing closing brace added once
    TEST_EQ(do_format<CharType>(loc, "End}}"), ascii_to<CharType>("End}"));
    // ...and twice when another trailing brace is added
    TEST_EQ(do_format<CharType>(loc, "End}}}"), ascii_to<CharType>("End}}"));
    // Escaped braces
    TEST_EQ(do_format<CharType>(loc, "Unexpected {{ in file"), ascii_to<CharType>("Unexpected { in file"));
    TEST_EQ(do_format<CharType>(loc, "Unexpected {{ in {1}#{2}", "f", 7), ascii_to<CharType>("Unexpected { in f#7"));
    TEST_EQ(do_format<CharType>(loc, "Unexpected }} in file"), ascii_to<CharType>("Unexpected } in file"));
    TEST_EQ(do_format<CharType>(loc, "Unexpected }} in {1}#{2}", "f", 9), ascii_to<CharType>("Unexpected } in f#9"));

    // format with multiple types
    TEST_EQ(do_format<CharType>(loc, "{1} {2}", "hello", 2), ascii_to<CharType>("hello 2"));

    // format with locale & encoding
    {
#if BOOST_LOCALE_ICU_VERSION >= 400
        const auto expected = boost::locale::conv::utf_to_utf<CharType>("10,00\xC2\xA0€");
#else
        const auto expected = boost::locale::conv::utf_to_utf<CharType>("10,00 €"); // LCOV_EXCL_LINE
#endif
        TEST_EQ(do_format<CharType>(loc, "{1,cur,locale=de_DE.UTF-8}", 10), expected);
    }

#define TEST_FORMAT_CLS(fmt_string, value, expected_str) \
    test_format_class_impl<CharType>(fmt_string, value, expected_str, loc, __LINE__)

    // Test different types and modifiers
    TEST_FORMAT_CLS("{1}", 1200.1, "1200.1");
    TEST_FORMAT_CLS("Test {1,num}", 1200.1, "Test 1,200.1");
    TEST_FORMAT_CLS("{{1}} {1,number}", 3200.4, "{1} 3,200.4");
    // placeholder in escaped braces, see issue #194
    TEST_FORMAT_CLS("{{{1}}}", "num", "{num}");
    TEST_FORMAT_CLS("{{{1}}}", 1200.1, "{1200.1}");

    TEST_FORMAT_CLS("{1,num=sci,p=3}", 13.1, "1.310E1");
    TEST_FORMAT_CLS("{1,num=scientific,p=3}", 13.1, "1.310E1");
    TEST_FORMAT_CLS("{1,num=fix,p=3}", 13.1, "13.100");
    TEST_FORMAT_CLS("{1,num=fixed,p=3}", 13.1, "13.100");
    TEST_FORMAT_CLS("{1,num=hex}", 0x1234, "1234");
    TEST_FORMAT_CLS("{1,num=oct}", 42, "52");
    TEST_FORMAT_CLS("{1,<,w=3,num}", -1, "-1 ");
    TEST_FORMAT_CLS("{1,>,w=3,num}", 1, "  1");
    TEST_FORMAT_CLS("{per,1}", 0.1, "10%");
    TEST_FORMAT_CLS("{percent,1}", 0.1, "10%");
    TEST_FORMAT_CLS("{1,cur}", 1234, "$1,234.00");
    TEST_FORMAT_CLS("{1,currency}", 1234, "$1,234.00");
    if(charset == "UTF-8") {
#if BOOST_LOCALE_ICU_VERSION >= 400
        TEST_FORMAT_CLS("{1,cur,locale=de_DE}", 10, "10,00\xC2\xA0€");
#else
        TEST_FORMAT_CLS("{1,cur,locale=de_DE}", 10, "10,00 €");                     // LCOV_EXCL_LINE
#endif
    }
#if BOOST_LOCALE_ICU_VERSION >= 402
    TEST_FORMAT_CLS("{1,cur=nat}", 1234, "$1,234.00");
    TEST_FORMAT_CLS("{1,cur=national}", 1234, "$1,234.00");
    TEST_FORMAT_CLS("{1,cur=iso}", 1234, get_icu_currency_iso(1234));
#endif
    TEST_FORMAT_CLS("{1,spell}", 10, "ten");
    TEST_FORMAT_CLS("{1,spellout}", 10, "ten");
#if 402 <= BOOST_LOCALE_ICU_VERSION && BOOST_LOCALE_ICU_VERSION < 408
    if(charset == "UTF-8") {
        TEST_FORMAT_CLS("{1,ord}", 1, "1\xcb\xa2\xe1\xb5\x97");
        TEST_FORMAT_CLS("{1,ordinal}", 1, "1\xcb\xa2\xe1\xb5\x97");
    }
#else
    TEST_FORMAT_CLS("{1,ord}", 1, "1st");
    TEST_FORMAT_CLS("{1,ordinal}", 1, "1st");
#endif

    // formatted time
    {
        boost::locale::time_zone::global("GMT+4:00");
        const time_t now = time(nullptr);
        char local_time_str[256], local_time_str_gmt2[256];
        time_t local_now = now + 3600 * 4;
        strftime(local_time_str, sizeof(local_time_str), "'%H:%M:%S'", gmtime_wrap(&local_now));
        local_now = now + 3600 * 2;
        strftime(local_time_str_gmt2, sizeof(local_time_str_gmt2), "'%H:%M:%S'", gmtime_wrap(&local_now));

        TEST_FORMAT_CLS("{1,ftime='''%H:%M:%S'''}", now, local_time_str);
        TEST_FORMAT_CLS("{1,strftime='''%H:%M:%S'''}", now, local_time_str);
        // 'local' has no impact on str(), uses global timezone
        TEST_FORMAT_CLS("{1,local,ftime='''%H:%M:%S'''}", now, local_time_str);

        std::basic_ostringstream<CharType> ss;
        ss.imbue(loc);
        ss << as::time_zone("GMT+02:00");
        format_type fmt_stream(ascii_to<CharType>("{1,ftime='''%H:%M:%S'''}"));      // Use timezone of stream
        format_type fmt_local(ascii_to<CharType>("{1,local,ftime='''%H:%M:%S'''}")); // Use global timezone
        empty_stream(ss) << fmt_stream % now;
        TEST_EQ(ss.str(), to<CharType>(local_time_str_gmt2));
        empty_stream(ss) << fmt_local % now;
        TEST_EQ(ss.str(), to<CharType>(local_time_str));
    }

    time_t a_date = 3600 * 24 * (31 + 4); // Feb 5th
    time_t a_time = 3600 * 15 + 60 * 33;  // 15:33:05
    time_t a_timesec = 13;
    time_t a_datetime = a_date + a_time + a_timesec;
    const std::string icu_time_def = get_ICU_time(as::time, a_datetime);
    const std::string icu_time_short = get_ICU_time(as::time_short, a_datetime);
    const std::string icu_time_medium = get_ICU_time(as::time_medium, a_datetime);
    const std::string icu_time_long = get_ICU_time(as::time_long, a_datetime);
    const std::string icu_time_full = get_ICU_time(as::time_full, a_datetime);
    const std::string icu_date_short = get_ICU_date(as::date_short, a_datetime);
    const std::string icu_date_medium = get_ICU_date(as::date_medium, a_datetime);
    const std::string icu_date_long = get_ICU_date(as::date_long, a_datetime);
    const std::string icu_date_full = get_ICU_date(as::date_full, a_datetime);
    const std::string icu_datetime_def = get_ICU_datetime(as::time, a_datetime);
    const std::string icu_datetime_short = get_ICU_datetime(as::time_short, a_datetime);
    const std::string icu_datetime_medium = get_ICU_datetime(as::time_medium, a_datetime);
    const std::string icu_datetime_long = get_ICU_datetime(as::time_long, a_datetime);
    const std::string icu_datetime_full = get_ICU_datetime(as::time_full, a_datetime);
    // Sanity check
    TEST_EQ(icu_date_full, "Thursday, February 5, 1970");

    TEST_FORMAT_CLS("{1,date,gmt}", a_datetime, "Feb 5, 1970");
    TEST_FORMAT_CLS("{1,time,gmt}", a_datetime, icu_time_def);
    TEST_FORMAT_CLS("{1,datetime,gmt}", a_datetime, icu_datetime_def);
    TEST_FORMAT_CLS("{1,dt,gmt}", a_datetime, icu_datetime_def);
    // With length modifier
    TEST_FORMAT_CLS("{1,time=short,gmt}", a_datetime, icu_time_short);
    TEST_FORMAT_CLS("{1,time=s,gmt}", a_datetime, icu_time_short);
    TEST_FORMAT_CLS("{1,time=medium,gmt}", a_datetime, icu_time_medium);
    TEST_FORMAT_CLS("{1,time=m,gmt}", a_datetime, icu_time_medium);
    TEST_FORMAT_CLS("{1,time=long,gmt}", a_datetime, icu_time_long);
    TEST_FORMAT_CLS("{1,time=l,gmt}", a_datetime, icu_time_long);
    TEST_FORMAT_CLS("{1,time=full,gmt}", a_datetime, icu_time_full);
    TEST_FORMAT_CLS("{1,time=f,gmt}", a_datetime, icu_time_full);
    TEST_FORMAT_CLS("{1,date=short,gmt}", a_datetime, icu_date_short);
    TEST_FORMAT_CLS("{1,date=s,gmt}", a_datetime, icu_date_short);
    TEST_FORMAT_CLS("{1,date=medium,gmt}", a_datetime, icu_date_medium);
    TEST_FORMAT_CLS("{1,date=m,gmt}", a_datetime, icu_date_medium);
    TEST_FORMAT_CLS("{1,date=long,gmt}", a_datetime, icu_date_long);
    TEST_FORMAT_CLS("{1,date=l,gmt}", a_datetime, icu_date_long);
    TEST_FORMAT_CLS("{1,date=full,gmt}", a_datetime, icu_date_full);
    TEST_FORMAT_CLS("{1,date=f,gmt}", a_datetime, icu_date_full);
    // Handle timezones and reuse of arguments
    const std::string icu_time_short2 = get_ICU_time(as::time_short, a_datetime, "GMT+01:00");
    TEST_FORMAT_CLS("{1,time=s,gmt};{1,time=s,timezone=GMT+01:00}", a_datetime, icu_time_short + ";" + icu_time_short2);
    TEST_FORMAT_CLS("{1,time=s,gmt};{1,time=s,tz=GMT+01:00}", a_datetime, icu_time_short + ";" + icu_time_short2);
    // Handle single quotes
    TEST_FORMAT_CLS("{1,gmt,ftime='%H'''}", a_datetime, "15'");
    TEST_FORMAT_CLS("{1,gmt,ftime='''%H'}", a_datetime, "'15");
    TEST_FORMAT_CLS("{1,gmt,ftime='%H o''clock'}", a_datetime, "15 o'clock");

    // Test not a year of the week
    a_datetime = 1388491200; // 2013-12-31 12:00 - check we don't use week of year

    TEST_FORMAT_CLS("{1,gmt,ftime='%Y'}", a_datetime, "2013");
    TEST_FORMAT_CLS("{1,gmt,ftime='%y'}", a_datetime, "13");
    TEST_FORMAT_CLS("{1,gmt,ftime='%D'}", a_datetime, "12/31/13");
}

BOOST_LOCALE_DISABLE_UNREACHABLE_CODE_WARNING
void test_main(int argc, char** argv)
{
    if(argc == 2)
        message_path = argv[1];

#ifndef BOOST_LOCALE_WITH_ICU
    std::cout << "ICU is not build... Skipping\n";
    return;
#endif
    boost::locale::time_zone::global("GMT+4:00");
    std::cout << "Testing char, UTF-8" << std::endl;
    test_manip<char>();
    test_format_class<char>();
    std::cout << "Testing char, ISO8859-1" << std::endl;
    test_manip<char>("ISO8859-1");
    test_format_class<char>("ISO8859-1");

    std::cout << "Testing wchar_t" << std::endl;
    test_manip<wchar_t>();
    test_format_class<wchar_t>();

#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
    std::cout << "Testing char16_t" << std::endl;
    test_manip<char16_t>();
    test_format_class<char16_t>();
#endif

#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
    std::cout << "Testing char32_t" << std::endl;
    test_manip<char32_t>();
    test_format_class<char32_t>();
#endif
}

// boostinspect:noascii
// boostinspect:nominmax
