//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_LOCALE_TEST_TOOLS_HPP
#define BOOST_LOCALE_TEST_TOOLS_HPP

#include <boost/locale/encoding.hpp>
#include "boostLocale/test/posix_tools.hpp"
#include "boostLocale/test/unit_test.hpp"
#include <cstdio>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#ifndef BOOST_LOCALE_NO_WINAPI_BACKEND
#    include "../src/boost/locale/win32/lcid.hpp"
#else
#    include <boost/core/ignore_unused.hpp>
#endif
#if BOOST_LOCALE_USE_WIN32_API
#    ifndef NOMINMAX
#        define NOMINMAX
#    endif
#    include <windows.h>
bool hasWinCodepage(unsigned codepage)
{
    return IsValidCodePage(codepage) != 0;
}
#else
bool hasWinCodepage(unsigned)
{
    return false;
}
#endif

#if defined(BOOST_MSVC) && BOOST_MSVC < 1700
#    pragma warning(disable : 4428) // universal-character-name encountered in source
#endif

class remove_file_on_exit {
    std::string filename_;

public:
    explicit remove_file_on_exit(const std::string& filename) : filename_(filename) {}
    ~remove_file_on_exit() { std::remove(filename_.c_str()); }
};

inline unsigned utf8_next(const std::string& s, unsigned& pos)
{
    unsigned c = static_cast<unsigned char>(s[pos++]);
    unsigned l;
    if(c <= 127)
        return c;
    else if(c <= 193)
        throw std::logic_error("Invalid UTF8"); // LCOV_EXCL_LINE
    else if(c <= 223)
        l = 1;
    else if(c <= 239)
        l = 2;
    else if(c <= 244)
        l = 3;
    else
        throw std::logic_error("Invalid UTF8"); // LCOV_EXCL_LINE

    c &= (1 << (6 - l)) - 1;

    switch(l) {
        case 3: c = (c << 6) | (static_cast<unsigned char>(s[pos++]) & 0x3F); BOOST_FALLTHROUGH;
        case 2: c = (c << 6) | (static_cast<unsigned char>(s[pos++]) & 0x3F); BOOST_FALLTHROUGH;
        case 1: c = (c << 6) | (static_cast<unsigned char>(s[pos++]) & 0x3F);
    }
    return c;
}

/// Convert an UTF encoded string to an UTF-8 encoded string
template<typename C>
std::string to_utf8(const std::basic_string<C>& utf_string)
{
    return boost::locale::conv::utf_to_utf<char>(utf_string);
}
std::string to_utf8(const std::string& utf_string)
{
    return utf_string;
}

/// Convert/decode an UTF-8 encoded string to the given char type
/// For `char` this will be Latin1, otherwise UTF-16/UTF-32
template<typename Char>
std::basic_string<Char> to(const std::string& utf8)
{
    std::basic_string<Char> out;
    for(unsigned i = 0; i < utf8.size();) {
        const unsigned prev = i;
        unsigned point = utf8_next(utf8, i);
        BOOST_LOCALE_START_CONST_CONDITION
        if(sizeof(Char) == 1 && point > 255) {
            std::ostringstream ss;
            ss << "Can't convert codepoint U" << std::hex << point << "("
               << std::string(utf8.begin() + prev, utf8.begin() + i) << ") to Latin1";
            throw std::logic_error(ss.str());
        } else if(sizeof(Char) == 2 && point > 0xFFFF) { // Deal with surrogates
            point -= 0x10000;
            out += static_cast<Char>(0xD800 | (point >> 10));
            out += static_cast<Char>(0xDC00 | (point & 0x3FF));
            continue;
        }
        BOOST_LOCALE_END_CONST_CONDITION
        out += static_cast<Char>(point);
    }
    return out;
}

#ifndef BOOST_LOCALE_NO_CXX20_STRING8
template<>
std::basic_string<char8_t> to(const std::string& utf8)
{
    return std::basic_string<char8_t>(utf8.begin(), utf8.end());
}
#endif

/// Convert an ASCII string to the given char type (i.e. copy only)
template<typename Char, size_t size>
inline std::basic_string<Char> ascii_to(const char (&str)[size])
{
    return std::basic_string<Char>(str, str + size - 1);
}

/// Convert an UTF-8 encoded string to another UTF encoding
/// or to a narrow string encoded using the given locale
template<typename Char>
std::basic_string<Char> to_correct_string(const std::string& utf8_str, std::locale /*l*/)
{
    return to<Char>(utf8_str);
}

/// Specialization to convert an UTF-8 encoded string to a locale specific encoded string
template<>
inline std::string to_correct_string(const std::string& utf8_str, std::locale l)
{
    return boost::locale::conv::from_utf(utf8_str, l);
}

bool has_std_locale(const char* name)
{
    try {
        std::locale tmp(name);
        return true;
    } catch(...) {
        return false;
    }
}

bool has_win_locale(const std::string& locale_name)
{
#ifdef BOOST_LOCALE_NO_WINAPI_BACKEND
    boost::ignore_unused(locale_name); // LCOV_EXCL_LINE
    return false;                      // LCOV_EXCL_LINE
#else
    return boost::locale::impl_win::locale_to_lcid(locale_name) != 0;
#endif
}

/// Clear a string stream and return it
template<class T>
T& empty_stream(T& s)
{
    s.str(std::basic_string<typename T::char_type>());
    s.clear();
    return s;
}

inline bool test_std_supports_SJIS_codecvt(const std::string& locale_name)
{
    const std::string file_path = boost::locale::test::exe_name + "-test-siftjis.txt";
    remove_file_on_exit _(file_path);
    {
        // Japan in Shift JIS/cp932
        const char* japan_932 = "\x93\xfa\x96\x7b";
        std::ofstream f(file_path, std::ios::binary);
        f << japan_932;
    }
    bool res = true;
    try {
        std::wfstream test;
        test.imbue(std::locale(locale_name));
        test.open(file_path);
        // Japan in Unicode
        const std::wstring cmp = L"\u65e5\u672c";
        std::wstring ref;
        res = (test >> ref) && (ref == cmp);
    } catch(const std::exception&) {
        res = false;
    }
    return res;
}

std::string get_std_name(const std::string& name, std::string* real_name = nullptr)
{
    if(has_std_locale(name.c_str())) {
        if(real_name)
            *real_name = name;
        return name;
    }

#if BOOST_LOCALE_USE_WIN32_API
    const bool utf8 = name.find("UTF-8") != std::string::npos;

    if(name == "en_US.UTF-8" || name == "en_US.ISO8859-1") {
        if(has_std_locale("English_United States.1252")) {
            if(real_name)
                *real_name = "English_United States.1252";
            return utf8 ? name : "en_US.windows-1252";
        }
        return "";
    } else if(name == "he_IL.UTF-8" || name == "he_IL.ISO8859-8") {
        if(has_std_locale("Hebrew_Israel.1255")) {
            if(real_name)
                *real_name = "Hebrew_Israel.1255";
            return utf8 ? name : "he_IL.windows-1255";
        }
    } else if(name == "ru_RU.UTF-8") {
        if(has_std_locale("Russian_Russia.1251")) {
            if(real_name)
                *real_name = "Russian_Russia.1251";
            return name;
        }
    } else if(name == "tr_TR.UTF-8") {
        if(has_std_locale("Turkish_Turkey.1254")) {
            if(real_name)
                *real_name = "Turkish_Turkey.1254";
            return name;
        }
    }
    if(name == "ja_JP.SJIS") {
        if(has_std_locale("Japanese_Japan.932")) {
            if(real_name)
                *real_name = "Japanese_Japan.932";
            return name;
        }
        return "";
    }
#endif
    return "";
}

char* make2(unsigned v)
{
    static unsigned char buf[3] = {0};
    buf[0] = static_cast<unsigned char>(0xC0 | (v >> 6));
    buf[1] = static_cast<unsigned char>(0x80 | (v & 0x3F));
    return reinterpret_cast<char*>(buf);
}

char* make3(unsigned v)
{
    static unsigned char buf[4] = {0};
    buf[0] = static_cast<unsigned char>(0xE0 | ((v >> 12)));
    buf[1] = static_cast<unsigned char>(0x80 | ((v >> 6) & 0x3F));
    buf[2] = static_cast<unsigned char>(0x80 | ((v >> 0) & 0x3F));
    return reinterpret_cast<char*>(buf);
}

char* make4(unsigned v)
{
    static unsigned char buf[5] = {0};
    buf[0] = static_cast<unsigned char>(0xF0 | ((v >> 18)));
    buf[1] = static_cast<unsigned char>(0x80 | ((v >> 12) & 0x3F));
    buf[2] = static_cast<unsigned char>(0x80 | ((v >> 6) & 0x3F));
    buf[3] = static_cast<unsigned char>(0x80 | ((v >> 0) & 0x3F));
    return reinterpret_cast<char*>(buf);
}

#ifdef _MSC_VER
#    pragma warning(push)
#    pragma warning(disable : 4996) //"This function or variable may be unsafe"
#endif
#if defined(__clang__)
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif
/// Wrapper for std::gmtime avoiding warning 4996 on MSVC/clang-cl:
inline std::tm* gmtime_wrap(const std::time_t* time)
{
    return std::gmtime(time);
}
/// Wrapper for std::localtime avoiding warning 4996 on MSVC/clang-cl
inline std::tm* localtime_wrap(const std::time_t* time)
{
    return std::localtime(time);
}
#if defined(__clang__)
#    pragma clang diagnostic pop
#endif
#ifdef _MSC_VER
#    pragma warning(pop)
#endif

#endif
