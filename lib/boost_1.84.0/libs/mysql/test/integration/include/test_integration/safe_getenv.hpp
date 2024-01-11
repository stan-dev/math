//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_SAFE_GETENV_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_SAFE_GETENV_HPP

#include <cstdlib>
#include <string>

namespace boost {
namespace mysql {
namespace test {

inline std::string safe_getenv(const char* name, const char* default_value)
{
    // MSVC doesn't like getenv
#ifdef BOOST_MSVC
#pragma warning(push)
#pragma warning(disable : 4996)
#endif
    const char* res = std::getenv(name);
#ifdef BOOST_MSVC
#pragma warning(pop)
#endif
    return res ? res : default_value;
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
