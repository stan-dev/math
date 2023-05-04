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

#ifdef BOOST_LOCALE_WITH_ICU
std::string get_icu_gmt_name(icu::TimeZone::EDisplayType style)
{
    icu::UnicodeString tmp;
    return from_icu_string(icu::TimeZone::getGMT()->getDisplayName(false, style, get_icu_test_locale(), tmp));
}

// This changes between ICU versions, e.g. "GMT" or "Greenwich Mean Time"
const std::string icu_full_gmt_name = get_icu_gmt_name(icu::TimeZone::EDisplayType::LONG);
// e.g. "GMT", "GMT+00:00"
const std::string icu_gmt_name = get_icu_gmt_name(icu::TimeZone::EDisplayType::SHORT);
#else
const std::string icu_full_gmt_name, icu_gmt_name;
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

#define TEST_MIN_MAX_FMT(type, minval, maxval)                      \
    TEST_FMT(as::number, std::numeric_limits<type>::min(), minval); \
    TEST_FMT(as::number, std::numeric_limits<type>::max(), maxval)

#define TEST_MIN_MAX_PARSE(type, minval, maxval)                      \
    TEST_PARSE(as::number, minval, std::numeric_limits<type>::min()); \
    TEST_PARSE(as::number, maxval, std::numeric_limits<type>::max())

#define TEST_MIN_MAX(type, minval, maxval)  \
    TEST_MIN_MAX_FMT(type, minval, maxval); \
    TEST_MIN_MAX_PARSE(type, minval, maxval)

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

template<typename CharType>
void test_manip(std::string e_charset = "UTF-8")
{
    using string_type = std::basic_string<CharType>;
    boost::locale::generator g;
    std::locale loc = g(test_locale_name + "." + e_charset);

    TEST_FMT_PARSE_1(as::posix, 1200.1, "1200.1");
    TEST_FMT_PARSE_1(as::number, 1200.1, "1,200.1");
    TEST_FMT(as::number << std::setfill(CharType('_')) << std::setw(6), 1534, "_1,534");
    TEST_FMT(as::number << std::left << std::setfill(CharType('_')) << std::setw(6), 1534, "1,534_");

    // Ranges
    TEST_MIN_MAX(int16_t, "-32,768", "32,767");
    TEST_MIN_MAX(uint16_t, "0", "65,535");
    TEST_PARSE_FAILS(as::number, "-1", uint16_t);
    if(stdlib_correctly_errors_on_out_of_range_int16()) {
        TEST_PARSE_FAILS(as::number, "65,535", int16_t);
    }

    TEST_MIN_MAX(int32_t, "-2,147,483,648", "2,147,483,647");
    TEST_MIN_MAX(uint32_t, "0", "4,294,967,295");
    TEST_PARSE_FAILS(as::number, "-1", uint32_t);
    TEST_PARSE_FAILS(as::number, "4,294,967,295", int32_t);

    TEST_MIN_MAX(int64_t, "-9,223,372,036,854,775,808", "9,223,372,036,854,775,807");
    // ICU does not support uint64, but we have a fallback to format it at least
    TEST_MIN_MAX_FMT(uint64_t, "0", "18446744073709551615");
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
    if(e_charset == "UTF-8") {
        TEST_FMT(as::ordinal, 1, "1\xcb\xa2\xe1\xb5\x97"); // 1st with st as ligatures
    }
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

    TEST_FMT_PARSE_2_2(as::time, as::gmt, a_datetime, "3:33:13 PM", a_time + a_timesec);
    TEST_FMT_PARSE_3_2(as::time, as::time_short, as::gmt, a_datetime, "3:33 PM", a_time);
    TEST_FMT_PARSE_3_2(as::time, as::time_medium, as::gmt, a_datetime, "3:33:13 PM", a_time + a_timesec);
    TEST_FMT_PARSE_3_2(as::time, as::time_long, as::gmt, a_datetime, "3:33:13 PM " + icu_gmt_name, a_time + a_timesec);
    // ICU 4.8.0 has a bug which makes parsing the full time fail when anything follows the time zone
#if BOOST_LOCALE_ICU_VERSION_EXACT != 40800
    TEST_FMT_PARSE_3_2(as::time,
                       as::time_full,
                       as::gmt,
                       a_datetime,
                       "3:33:13 PM " + icu_full_gmt_name,
                       a_time + a_timesec);
#endif
    TEST_PARSE_FAILS(as::time, "AM", double);

    TEST_FMT_PARSE_2_2(as::time, as::time_zone("GMT+01:00"), a_datetime, "4:33:13 PM", a_time + a_timesec);
    TEST_FMT_PARSE_3_2(as::time, as::time_short, as::time_zone("GMT+01:00"), a_datetime, "4:33 PM", a_time);
    TEST_FMT_PARSE_3_2(as::time,
                       as::time_medium,
                       as::time_zone("GMT+01:00"),
                       a_datetime,
                       "4:33:13 PM",
                       a_time + a_timesec);

#if U_ICU_VERSION_MAJOR_NUM >= 51
#    define GMT_P100 "GMT+1"
#else
#    define GMT_P100 "GMT+01:00"
#endif

#if U_ICU_VERSION_MAJOR_NUM >= 50
#    define ICU_COMMA ","
#    define ICUAT " at"
#else
#    define ICU_COMMA ""
#    define ICUAT ""
#endif

    TEST_FMT_PARSE_3_2(as::time,
                       as::time_long,
                       as::time_zone("GMT+01:00"),
                       a_datetime,
                       "4:33:13 PM " GMT_P100,
                       a_time + a_timesec);
#if !(BOOST_LOCALE_ICU_VERSION == 308 && defined(__CYGWIN__)) // Known failure due to ICU issue
    TEST_FMT_PARSE_3_2(as::time,
                       as::time_full,
                       as::time_zone("GMT+01:00"),
                       a_datetime,
                       "4:33:13 PM GMT+01:00",
                       a_time + a_timesec);
#endif

    TEST_FMT_PARSE_2(as::datetime, as::gmt, a_datetime, "Feb 5, 1970" ICU_COMMA " 3:33:13 PM");
    TEST_FMT_PARSE_4_2(as::datetime,
                       as::date_short,
                       as::time_short,
                       as::gmt,
                       a_datetime,
                       "2/5/70" ICU_COMMA " 3:33 PM",
                       a_date + a_time);
    TEST_FMT_PARSE_4(as::datetime,
                     as::date_medium,
                     as::time_medium,
                     as::gmt,
                     a_datetime,
                     "Feb 5, 1970" ICU_COMMA " 3:33:13 PM");
    TEST_FMT_PARSE_4(as::datetime,
                     as::date_long,
                     as::time_long,
                     as::gmt,
                     a_datetime,
                     "February 5, 1970" ICUAT " 3:33:13 PM " + icu_gmt_name);
#if BOOST_LOCALE_ICU_VERSION_EXACT != 40800
    // ICU 4.8.0 has a bug which makes parsing the full time fail when anything follows the time zone
    TEST_FMT_PARSE_4(as::datetime,
                     as::date_full,
                     as::time_full,
                     as::gmt,
                     a_datetime,
                     "Thursday, February 5, 1970" ICUAT " 3:33:13 PM " + icu_full_gmt_name);
#endif

    const std::pair<char, std::string> mark_test_cases[] = {
      std::make_pair('a', "Thu"),
      std::make_pair('A', "Thursday"),
      std::make_pair('b', "Feb"),
      std::make_pair('B', "February"),
      std::make_pair('c', "Thursday, February 5, 1970" ICUAT " 3:33:13 PM " + icu_full_gmt_name),
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
      std::make_pair('X', "3:33:13 PM"),
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
        TEST_EQ(ss.str(), to<CharType>(mark_result.second));
    }

    {
        const time_t now = time(0);
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
        ss.str(string_type()); // Clear
        // We can manually tell it to use the local time zone
        ss << as::ftime(format_string) << as::local_time << now;
        TEST_EQ(ss.str(), to<CharType>(time_str));
        ss.str(string_type()); // Clear
        // Or e.g. GMT
        ss << as::ftime(format_string) << as::gmt << now;
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
            ss << format_type(fmt_string) % 1 % 2 % 3;
            TEST_EQ(ss.str(), expected);

            // Stream translated output
            ss.str(string_type());
            ss << format_type(boost::locale::translate(fmt_string)) % 1 % 2 % 3;
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
            int i1 = 1, i2 = 2, i3 = 3;
            fmt2 % i1 % i2 % i3;
            fmt2 = format_type(ascii_to<CharType>("{1}"));
            TEST_EQ(fmt2.str(), ascii_to<CharType>("")); // No bound value
            TEST_EQ((fmt2 % 42).str(), ascii_to<CharType>("42"));

            fmt2 = format_type(hello);
            TEST_EQ(fmt2.str(), hello); // Not translated
            fmt2 = format_type(boost::locale::translate(hello));
            TEST_EQ(fmt2.str(), hello_he); // Translated
        }
        // Restore
        std::locale::global(old_locale);
    }

    // Not passed placeholders are removed
    TEST_EQ((format_type(ascii_to<CharType>("{1}{3}{2}")) % "hello" % "world").str(loc),
            ascii_to<CharType>("helloworld"));
    TEST_EQ(format_type(ascii_to<CharType>("{1}")).str(loc), ascii_to<CharType>(""));
    // Unexpected closing brace and other chars are ignored
    TEST_EQ(format_type(ascii_to<CharType>(" = , } 3")).str(loc), ascii_to<CharType>(" = , } 3"));
    // Trailing opening brace is ignored
    TEST_EQ(format_type(ascii_to<CharType>("End {")).str(loc), ascii_to<CharType>("End "));
    // Trailing closing brace is added like any other char
    TEST_EQ(format_type(ascii_to<CharType>("End}")).str(loc), ascii_to<CharType>("End}"));
    // Escaped trailing closing brace added once
    TEST_EQ(format_type(ascii_to<CharType>("End}}")).str(loc), ascii_to<CharType>("End}"));
    // ...and twice when another trailing brace is added
    TEST_EQ(format_type(ascii_to<CharType>("End}}}")).str(loc), ascii_to<CharType>("End}}"));

    // format with multiple types
    TEST_EQ((format_type(ascii_to<CharType>("{1} {2}")) % "hello" % 2).str(loc), ascii_to<CharType>("hello 2"));

#define TEST_FORMAT_CLS(fmt_string, value, expected_str) \
    test_format_class_impl<CharType>(fmt_string, value, expected_str, loc, __LINE__)

    // Test different types and modifiers
    TEST_FORMAT_CLS("{1}", 1200.1, "1200.1");
    TEST_FORMAT_CLS("Test {1,num}", 1200.1, "Test 1,200.1");
    TEST_FORMAT_CLS("{{}} {1,number}", 1200.1, "{} 1,200.1");
    TEST_FORMAT_CLS("{1,num=sci,p=3}", 13.1, "1.310E1");
    TEST_FORMAT_CLS("{1,num=scientific,p=3}", 13.1, "1.310E1");
    TEST_FORMAT_CLS("{1,num=fix,p=3}", 13.1, "13.100");
    TEST_FORMAT_CLS("{1,num=fixed,p=3}", 13.1, "13.100");
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
        TEST_FORMAT_CLS("{1,cur,locale=de_DE}", 10, "10,00 €");
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
        time_t now = time(0);
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
        ss << fmt_stream % now;
        TEST_EQ(ss.str(), to<CharType>(local_time_str_gmt2));
        ss.str(string_type());
        ss << fmt_local % now;
        TEST_EQ(ss.str(), to<CharType>(local_time_str));
    }

    time_t a_date = 3600 * 24 * (31 + 4); // Feb 5th
    time_t a_time = 3600 * 15 + 60 * 33;  // 15:33:05
    time_t a_timesec = 13;
    time_t a_datetime = a_date + a_time + a_timesec;
    TEST_FORMAT_CLS("{1,date,gmt}", a_datetime, "Feb 5, 1970");
    TEST_FORMAT_CLS("{1,time,gmt}", a_datetime, "3:33:13 PM");
    TEST_FORMAT_CLS("{1,datetime,gmt}", a_datetime, "Feb 5, 1970" ICU_COMMA " 3:33:13 PM");
    TEST_FORMAT_CLS("{1,dt,gmt}", a_datetime, "Feb 5, 1970" ICU_COMMA " 3:33:13 PM");
    // With length modifier
    TEST_FORMAT_CLS("{1,time=short,gmt}", a_datetime, "3:33 PM");
    TEST_FORMAT_CLS("{1,time=s,gmt}", a_datetime, "3:33 PM");
    TEST_FORMAT_CLS("{1,time=medium,gmt}", a_datetime, "3:33:13 PM");
    TEST_FORMAT_CLS("{1,time=m,gmt}", a_datetime, "3:33:13 PM");
    TEST_FORMAT_CLS("{1,time=long,gmt}", a_datetime, "3:33:13 PM " + icu_gmt_name);
    TEST_FORMAT_CLS("{1,time=l,gmt}", a_datetime, "3:33:13 PM " + icu_gmt_name);
    TEST_FORMAT_CLS("{1,date=full,gmt}", a_datetime, "Thursday, February 5, 1970");
    TEST_FORMAT_CLS("{1,date=f,gmt}", a_datetime, "Thursday, February 5, 1970");
    // Handle timezones and reuse of arguments
    TEST_FORMAT_CLS("{1,time=s,gmt};{1,time=s,timezone=GMT+01:00}", a_datetime, "3:33 PM;4:33 PM");
    TEST_FORMAT_CLS("{1,time=s,gmt};{1,time=s,tz=GMT+01:00}", a_datetime, "3:33 PM;4:33 PM");
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
