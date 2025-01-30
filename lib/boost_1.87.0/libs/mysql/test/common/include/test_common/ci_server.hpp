//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_CI_SERVER_HPP
#define BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_CI_SERVER_HPP

#include <boost/config.hpp>

#include <cstdlib>
#include <string>

// Constant and utilities to interact with the CI database server

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

// Interface
inline std::string get_hostname() { return safe_getenv("BOOST_MYSQL_SERVER_HOST", "127.0.0.1"); }
constexpr const char* integ_user = "integ_user";
constexpr const char* integ_passwd = "integ_password";
constexpr const char* integ_db = "boost_mysql_integtests";
constexpr const char* default_unix_path = "/var/run/mysqld/mysqld.sock";

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
