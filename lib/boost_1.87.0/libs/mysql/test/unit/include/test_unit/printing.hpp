//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_PRINTING_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_PRINTING_HPP

#include <iosfwd>

namespace boost {
namespace mysql {

// address_type
enum class address_type;
std::ostream& operator<<(std::ostream& os, address_type value);

namespace detail {

// capabilities
class capabilities;
std::ostream& operator<<(std::ostream& os, const capabilities& caps);

// db_flavor
enum class db_flavor;
std::ostream& operator<<(std::ostream& os, db_flavor value);

// resultset_encoding
enum class resultset_encoding;
std::ostream& operator<<(std::ostream& os, resultset_encoding t);

// results_iterator
class results_iterator;
std::ostream& operator<<(std::ostream& os, const results_iterator& it);

// next_action_type
enum class next_action_type;
std::ostream& operator<<(std::ostream& os, next_action_type t);

// pipeline_stage_kind
enum class pipeline_stage_kind;
std::ostream& operator<<(std::ostream& os, pipeline_stage_kind v);

// pipeline_request_stage
struct pipeline_request_stage;
bool operator==(const pipeline_request_stage& lhs, const pipeline_request_stage& rhs);
std::ostream& operator<<(std::ostream& os, pipeline_request_stage v);

// connection_status (pool)
enum class connection_status;
std::ostream& operator<<(std::ostream& os, connection_status v);

// collection_state (pool)
enum class collection_state;
std::ostream& operator<<(std::ostream& os, collection_state v);

// next_connection_action (pool)
enum class next_connection_action;
std::ostream& operator<<(std::ostream& os, next_connection_action v);

}  // namespace detail
}  // namespace mysql
}  // namespace boost

#endif
