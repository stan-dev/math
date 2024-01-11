//
// Copyright (c) 2015 Artyom Beilis (Tonkikh)
// Copyright (c) 2021 Alexander Grund
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#if defined(__GNUC__) && __GNUC__ >= 7
#pragma GCC diagnostic ignored "-Wattributes"
#endif
#include <boost/nowide/filesystem.hpp>

#include <boost/nowide/convert.hpp>
#include <boost/nowide/cstdio.hpp>
#include <boost/nowide/fstream.hpp>
#include <boost/nowide/utf/convert.hpp>
#include "test.hpp"

#include <iomanip> // Required for feature macro check below
// Conditional include to avoid warning/message
#if defined(__cpp_lib_quoted_string_io) && __cpp_lib_quoted_string_io >= 201304
#include <boost/nowide/quoted.hpp>
#endif

#include <sstream>
#include <type_traits>
#if defined(_MSC_VER)
#pragma warning(disable : 4714) // function marked as __forceinline not inlined
#endif
#include <boost/filesystem.hpp>

// Exclude apple as support there is target level specific -.-
#if defined(__cpp_lib_filesystem) && !defined(__APPLE__)
#include <filesystem>
#define BOOST_NOWIDE_TEST_STD_PATH
#endif
#if defined(__cpp_lib_experimental_filesystem)
#ifndef _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#endif
#include <experimental/filesystem>
#define BOOST_NOWIDE_TEST_STD_EXPERIMENTAL_PATH
#endif

template<typename T, typename = void>
struct is_istreamable : std::false_type
{};
using boost::nowide::detail::void_t;
template<typename T>
struct is_istreamable<T, void_t<decltype(std::declval<std::istream&>() >> std::declval<T>())>> : std::true_type
{};

template<typename T_Char>
std::string maybe_narrow(const std::basic_string<T_Char>& s)
{
    return boost::nowide::narrow(s);
}

const std::string& maybe_narrow(const std::string& s)
{
    return s;
}

template<class Path>
void test_fs_path_io(std::string utf8_name)
{
#if defined(__cpp_lib_quoted_string_io) && __cpp_lib_quoted_string_io >= 201304
    Path path(boost::nowide::utf::convert_string<typename Path::value_type>(utf8_name));
    // Get native and UTF-8/narrow name here as the Path ctor may change the string (e.g. slash substitution)
    const auto nativeName = path.native();
    utf8_name = maybe_narrow(nativeName);
    // Output
    std::ostringstream s, sRef;
    sRef << std::quoted(utf8_name);
    s << boost::nowide::quoted(path);
    TEST_EQ(s.str(), sRef.str());
    // const
    const Path constPath(path);
    s.str("");
    s << boost::nowide::quoted(constPath);
    TEST_EQ(s.str(), sRef.str());
    // Rvalue
    s.str("");
    s << boost::nowide::quoted(Path(path));
    TEST_EQ(s.str(), sRef.str());

    // Input
    std::istringstream sIn(sRef.str());
    Path pathOut;
    static_assert(is_istreamable<decltype(boost::nowide::quoted(pathOut))>::value, "!");
    sIn >> boost::nowide::quoted(pathOut);
    TEST_EQ(pathOut.native(), nativeName);
    // Can't read into a const path
    static_assert(!is_istreamable<decltype(boost::nowide::quoted(constPath))>::value, "!");
    // or an Rvalue
    static_assert(!is_istreamable<decltype(boost::nowide::quoted(Path(path)))>::value, "!");

    // Wide stream
    std::wostringstream ws, wsRef;
    wsRef << std::quoted(boost::nowide::widen(utf8_name));
    ws << boost::nowide::quoted(path);
    TEST_EQ(ws.str(), wsRef.str());
    std::wistringstream wsIn(wsRef.str());
    pathOut.clear();
    wsIn >> boost::nowide::quoted(pathOut);
    TEST_EQ(maybe_narrow(pathOut.native()), utf8_name);
#else
    (void)utf8_name; // Suppress unused warning
    std::cout << "Skipping tests for boost::nowide::quoted" << std::endl;
#endif
}

// coverity[root_function]
void test_main(int, char** argv, char**)
{
    boost::nowide::nowide_filesystem();
    const std::string prefix = argv[0];
    const std::string utf8_name =
      prefix + "\xf0\x9d\x92\x9e-\xD0\xBF\xD1\x80\xD0\xB8\xD0\xB2\xD0\xB5\xD1\x82-\xE3\x82\x84\xE3\x81\x82.txt";

    {
        boost::nowide::ofstream f(utf8_name.c_str());
        TEST(f);
        f << "Test" << std::endl;
    }

    TEST(boost::filesystem::is_regular_file(boost::nowide::widen(utf8_name)));
    TEST(boost::filesystem::is_regular_file(utf8_name));

    TEST(boost::nowide::remove(utf8_name.c_str()) == 0);

    TEST(!boost::filesystem::is_regular_file(boost::nowide::widen(utf8_name)));
    TEST(!boost::filesystem::is_regular_file(utf8_name));

    const boost::filesystem::path path = utf8_name;
    {
        boost::nowide::ofstream f(path);
        TEST(f);
        f << "Test" << std::endl;
        TEST(is_regular_file(path));
    }
    {
        boost::nowide::ifstream f(path);
        TEST(f);
        std::string test;
        f >> test;
        TEST(test == "Test");
    }
    {
        boost::nowide::fstream f(path);
        TEST(f);
        std::string test;
        f >> test;
        TEST(test == "Test");
    }
    boost::filesystem::remove(path);

    std::cout << "Testing boost::filesystem::path" << std::endl;
    test_fs_path_io<boost::filesystem::path>(utf8_name);
#ifdef BOOST_NOWIDE_TEST_STD_EXPERIMENTAL_PATH
    std::cout << "Testing std::experimental::filesystem::path" << std::endl;
    test_fs_path_io<std::experimental::filesystem::path>(utf8_name);
#endif
#ifdef BOOST_NOWIDE_TEST_STD_PATH
    std::cout << "Testing std::filesystem::path" << std::endl;
    test_fs_path_io<std::filesystem::path>(utf8_name);
#endif
}
