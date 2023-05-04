//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_UTILS_INCLUDE_ER_CONNECTION_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_UTILS_INCLUDE_ER_CONNECTION_HPP

#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/handshake_params.hpp>
#include <boost/mysql/metadata_mode.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/rows_view.hpp>
#include <boost/mysql/statement.hpp>

#include <forward_list>
#include <memory>

#include "network_result.hpp"

namespace boost {
namespace mysql {
namespace test {

using fv_list_it = std::forward_list<field_view>::const_iterator;

class er_network_variant;

class er_connection
{
public:
    virtual ~er_connection() {}
    virtual bool uses_ssl() const = 0;
    virtual bool is_open() const = 0;
    virtual void set_metadata_mode(metadata_mode v) = 0;
    virtual void physical_connect() = 0;     // used by fixture setup functions
    virtual void sync_close() noexcept = 0;  // used by fixture cleanup functions
    virtual er_network_variant& variant() const = 0;

    virtual network_result<void> connect(const handshake_params&) = 0;
    virtual network_result<void> handshake(const handshake_params&) = 0;
    virtual network_result<void> query(string_view query, results& result) = 0;
    virtual network_result<void> start_query(string_view query, execution_state& result) = 0;
    virtual network_result<statement> prepare_statement(string_view statement) = 0;
    virtual network_result<void> execute_statement(
        const statement& stmt,
        field_view fv1,
        field_view fv2,
        results& result
    ) = 0;
    virtual network_result<void> start_statement_execution(
        const statement& stmt,
        field_view fv1,
        field_view fv2,
        execution_state& st
    ) = 0;
    virtual network_result<void> start_statement_execution(
        const statement& stmt,
        fv_list_it params_first,
        fv_list_it params_last,
        execution_state& st
    ) = 0;
    virtual network_result<void> close_statement(statement&) = 0;
    virtual network_result<rows_view> read_some_rows(execution_state& st) = 0;
    virtual network_result<void> ping() = 0;
    virtual network_result<void> quit() = 0;
    virtual network_result<void> close() = 0;
};

using er_connection_ptr = std::unique_ptr<er_connection>;

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
