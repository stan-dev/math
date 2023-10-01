//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_ER_CONNECTION_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_ER_CONNECTION_HPP

#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/handshake_params.hpp>
#include <boost/mysql/metadata_mode.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/rows_view.hpp>
#include <boost/mysql/statement.hpp>
#include <boost/mysql/static_execution_state.hpp>
#include <boost/mysql/static_results.hpp>

#include <forward_list>
#include <memory>

#include "test_common/network_result.hpp"
#include "test_integration/static_rows.hpp"

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
    virtual network_result<void> execute(string_view, results&) = 0;
    virtual network_result<void> execute(bound_statement_tuple<std::tuple<field_view, field_view>>, results&) = 0;
    virtual network_result<void> execute(bound_statement_iterator_range<fv_list_it>, results&) = 0;
    virtual network_result<void> start_execution(string_view, execution_state&) = 0;
    virtual network_result<void> start_execution(bound_statement_tuple<std::tuple<field_view, field_view>>, execution_state&) = 0;
    virtual network_result<void> start_execution(bound_statement_iterator_range<fv_list_it>, execution_state&) = 0;
    virtual network_result<void> close_statement(statement&) = 0;
    virtual network_result<void> read_resultset_head(execution_state& st) = 0;
    virtual network_result<rows_view> read_some_rows(execution_state& st) = 0;
    virtual network_result<void> ping() = 0;
    virtual network_result<void> quit() = 0;
    virtual network_result<void> close() = 0;

#ifdef BOOST_MYSQL_CXX14
    using static_results_t = static_results<row_multifield, row_2fields, empty>;
    using static_state_t = static_execution_state<row_multifield, row_2fields, empty>;
    virtual network_result<void> execute(string_view, static_results_t&) = 0;
    virtual network_result<void> start_execution(string_view, static_state_t&) = 0;
    virtual network_result<void> read_resultset_head(static_state_t& st) = 0;
    virtual network_result<std::size_t> read_some_rows(static_state_t& st, boost::span<row_multifield>) = 0;
    virtual network_result<std::size_t> read_some_rows(static_state_t& st, boost::span<row_2fields>) = 0;
#endif
};

using er_connection_ptr = std::unique_ptr<er_connection>;

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
