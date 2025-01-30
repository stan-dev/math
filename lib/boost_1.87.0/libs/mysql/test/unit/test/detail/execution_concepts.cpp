//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/with_params.hpp>

#include <boost/mysql/detail/config.hpp>

#ifdef BOOST_MYSQL_HAS_CONCEPTS

#include <boost/mysql/statement.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/detail/execution_concepts.hpp>

#include <string>
#include <string_view>

using namespace boost::mysql;
using boost::mysql::detail::execution_state_type;
using boost::mysql::detail::is_execution_request;
using boost::mysql::detail::results_type;

namespace {

// -----
// is_execution_request
// -----

// strings
static_assert(is_execution_request<std::string>::value, "");
static_assert(is_execution_request<const std::string&>::value, "");
static_assert(is_execution_request<std::string&&>::value, "");
static_assert(is_execution_request<std::string&>::value, "");

static_assert(is_execution_request<string_view>::value, "");
static_assert(is_execution_request<const string_view&>::value, "");
static_assert(is_execution_request<string_view&&>::value, "");

static_assert(is_execution_request<std::string_view>::value, "");
static_assert(is_execution_request<const std::string_view&>::value, "");
static_assert(is_execution_request<std::string_view&&>::value, "");

static_assert(is_execution_request<const char*>::value, "");

static_assert(is_execution_request<const char[14]>::value, "");
static_assert(is_execution_request<const char (&)[14]>::value, "");

// tuple statements
using tup_type = std::tuple<int, const char*, std::nullptr_t>;
static_assert(is_execution_request<bound_statement_tuple<tup_type>>::value, "");
static_assert(is_execution_request<const bound_statement_tuple<tup_type>&>::value, "");
static_assert(is_execution_request<bound_statement_tuple<tup_type>&>::value, "");
static_assert(is_execution_request<bound_statement_tuple<tup_type>&&>::value, "");

static_assert(is_execution_request<bound_statement_tuple<std::tuple<>>>::value, "");

static_assert(!is_execution_request<tup_type>::value, "");
static_assert(!is_execution_request<const tup_type&>::value, "");
static_assert(!is_execution_request<tup_type&>::value, "");
static_assert(!is_execution_request<tup_type&&>::value, "");

// iterator range statements
static_assert(is_execution_request<bound_statement_iterator_range<field_view*>>::value, "");
static_assert(is_execution_request<const bound_statement_iterator_range<field_view*>&>::value, "");
static_assert(is_execution_request<bound_statement_iterator_range<field_view*>&>::value, "");
static_assert(is_execution_request<bound_statement_iterator_range<field_view*>&&>::value, "");

static_assert(!is_execution_request<field_view*>::value, "");

// with_params
static_assert(is_execution_request<with_params_t<>>::value, "");
static_assert(is_execution_request<with_params_t<int>>::value, "");
static_assert(is_execution_request<with_params_t<int, float>>::value, "");
static_assert(is_execution_request<with_params_t<const std::string&, float>>::value, "");

static_assert(is_execution_request<with_params_t<int>&>::value, "");
static_assert(is_execution_request<const with_params_t<int>&>::value, "");
static_assert(is_execution_request<with_params_t<int>&&>::value, "");
static_assert(is_execution_request<with_params_t<const std::string&>&&>::value, "");

// Other stuff
static_assert(!is_execution_request<field_view>::value, "");
static_assert(!is_execution_request<int>::value, "");
static_assert(!is_execution_request<std::nullptr_t>::value, "");

using row1 = std::tuple<int, float>;
using row2 = std::tuple<double>;

// -----
// execution_state_type
// -----
static_assert(!execution_state_type<std::string>, "");
static_assert(!execution_state_type<int>, "");
static_assert(!execution_state_type<results>, "");
static_assert(!execution_state_type<static_results<row1>>, "");
static_assert(!execution_state_type<static_results<row1, row2>>, "");
static_assert(execution_state_type<execution_state>, "");
static_assert(execution_state_type<static_execution_state<row1>>, "");
static_assert(execution_state_type<static_execution_state<row1, row2>>, "");

// -----
// results_type
// -----
static_assert(!results_type<std::string>, "");
static_assert(!results_type<int>, "");
static_assert(!results_type<execution_state>, "");
static_assert(!results_type<static_execution_state<row1>>, "");
static_assert(!results_type<static_execution_state<row1, row2>>, "");
static_assert(results_type<results>, "");
static_assert(results_type<static_results<row1>>, "");
static_assert(results_type<static_results<row1, row2>>, "");

}  // namespace

#endif
