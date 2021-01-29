//
//  Copyright (c) 2012 Artyom Beilis (Tonkikh)
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/nowide/convert.hpp>
#include <iostream>

#include "test.hpp"
#include "test_sets.hpp"

#if defined(BOOST_MSVC) && BOOST_MSVC < 1700
#pragma warning(disable : 4428) // universal-character-name encountered in source
#endif

std::wstring widen_buf_ptr(const std::string& s)
{
    wchar_t buf[50];
    TEST(boost::nowide::widen(buf, 50, s.c_str()) == buf);
    return buf;
}

std::string narrow_buf_ptr(const std::wstring& s)
{
    char buf[50];
    TEST(boost::nowide::narrow(buf, 50, s.c_str()) == buf);
    return buf;
}

std::wstring widen_buf_range(const std::string& s)
{
    wchar_t buf[50];
    TEST(boost::nowide::widen(buf, 50, s.c_str(), s.c_str() + s.size()) == buf);
    return buf;
}

std::string narrow_buf_range(const std::wstring& s)
{
    char buf[50];
    TEST(boost::nowide::narrow(buf, 50, s.c_str(), s.c_str() + s.size()) == buf);
    return buf;
}

std::wstring widen_raw_string(const std::string& s)
{
    return boost::nowide::widen(s.c_str());
}

std::string narrow_raw_string(const std::wstring& s)
{
    return boost::nowide::narrow(s.c_str());
}

std::wstring widen_raw_string_and_size(const std::string& s)
{
    // Remove NULL
    const std::string s2 = s + "DummyData";
    return boost::nowide::widen(s2.c_str(), s.size());
}

std::string narrow_raw_string_and_size(const std::wstring& s)
{
    // Remove NULL
    const std::wstring s2 = s + L"DummyData";
    return boost::nowide::narrow(s2.c_str(), s.size());
}

void test_main(int, char**, char**)
{
    std::string hello = "\xd7\xa9\xd7\x9c\xd7\x95\xd7\x9d";
    std::wstring whello = L"\u05e9\u05dc\u05d5\u05dd";
    std::wstring whello_3e = L"\u05e9\u05dc\u05d5\ufffd";
    std::wstring whello_3 = L"\u05e9\u05dc\u05d5";

    std::cout << "- boost::nowide::widen" << std::endl;
    {
        const char* b = hello.c_str();
        const char* e = b + hello.size();
        wchar_t buf[6] = {0, 0, 0, 0, 0, 1};
        TEST(boost::nowide::widen(buf, 5, b, e) == buf);
        TEST(buf == whello);
        TEST(buf[5] == 1);
        TEST(boost::nowide::widen(buf, 4, b, e) == 0);
        TEST(boost::nowide::widen(buf, 5, b, e - 1) == buf);
        TEST(buf == whello_3e);
        TEST(boost::nowide::widen(buf, 5, b, e - 2) == buf);
        TEST(buf == whello_3);
        TEST(boost::nowide::widen(buf, 5, b, b) == buf && buf[0] == 0);
        TEST(boost::nowide::widen(buf, 5, b, b + 2) == buf && buf[1] == 0 && buf[0] == whello[0]);
    }
    std::cout << "- boost::nowide::narrow" << std::endl;
    {
        const wchar_t* b = whello.c_str();
        const wchar_t* e = b + whello.size(); //-V594
        char buf[10] = {0};
        buf[9] = 1;
        TEST(boost::nowide::narrow(buf, 9, b, e) == buf);
        TEST(buf == hello);
        TEST(buf[9] == 1);
        TEST(boost::nowide::narrow(buf, 8, b, e) == 0);
        TEST(boost::nowide::narrow(buf, 7, b, e - 1) == buf);
        TEST(buf == hello.substr(0, 6));
    }

    std::cout << "- (output_buffer, buffer_size, input_raw_string)" << std::endl;
    run_all(widen_buf_ptr, narrow_buf_ptr);
    std::cout << "- (output_buffer, buffer_size, input_raw_string, string_len)" << std::endl;
    run_all(widen_buf_range, narrow_buf_range);
    std::cout << "- (input_raw_string)" << std::endl;
    run_all(widen_raw_string, narrow_raw_string);
    std::cout << "- (input_raw_string, size)" << std::endl;
    run_all(widen_raw_string_and_size, narrow_raw_string_and_size);
    std::cout << "- (const std::string&)" << std::endl;
    run_all(boost::nowide::widen, boost::nowide::narrow);
}
