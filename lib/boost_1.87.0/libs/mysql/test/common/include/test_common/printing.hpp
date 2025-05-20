//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_PRINTING_HPP
#define BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_PRINTING_HPP

#include <iosfwd>

namespace boost {
namespace mysql {

enum class client_errc;
std::ostream& operator<<(std::ostream& os, client_errc v);

enum class common_server_errc;
std::ostream& operator<<(std::ostream& os, common_server_errc v);

class diagnostics;
std::ostream& operator<<(std::ostream& os, const diagnostics& v);

class row_view;
std::ostream& operator<<(std::ostream& os, const row_view& v);

class row;
std::ostream& operator<<(std::ostream& os, const row& v);

enum class metadata_mode;
std::ostream& operator<<(std::ostream& os, metadata_mode v);

enum class ssl_mode;
std::ostream& operator<<(std::ostream& os, ssl_mode v);

struct character_set;
bool operator==(const character_set& lhs, const character_set& rhs);
std::ostream& operator<<(std::ostream& os, const character_set& v);

}  // namespace mysql
}  // namespace boost

#endif
