//
//  Copyright (c) 2020 Alexander Grund
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/nowide/detail/is_path.hpp>

#include "test.hpp"
#include <iostream>

#ifdef __has_include
#if __has_include(<version>)
#include <version>
#endif
#endif

// Exclude apple as support there is target level specific -.-
#if defined(__cpp_lib_filesystem) && !defined(__APPLE__)
#include <filesystem>
#define BOOST_NOWIDE_TEST_SFS_PATH
#endif

#ifdef BOOST_NOWIDE_TEST_BFS_PATH
#include <boost/filesystem/path.hpp>
#endif

void test_main(int, char**, char**)
{
#ifdef BOOST_NOWIDE_TEST_SFS_PATH
    std::cout << "Testing std::filesystem::path" << std::endl;
    static_assert(boost::nowide::detail::is_path<std::filesystem::path>::value, "!");
#endif
#ifdef BOOST_NOWIDE_TEST_BFS_PATH
    std::cout << "Testing boost::filesystem::path" << std::endl;
    static_assert(boost::nowide::detail::is_path<boost::filesystem::path>::value, "!");
#endif
}
