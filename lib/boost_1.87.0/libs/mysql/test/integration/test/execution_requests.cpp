//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/field.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/format_sql.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/sequence.hpp>
#include <boost/mysql/with_params.hpp>

#include <boost/asio/deferred.hpp>
#include <boost/optional/optional.hpp>
#include <boost/test/unit_test.hpp>

#include <forward_list>
#include <functional>
#include <string>
#include <vector>

#include "test_common/create_basic.hpp"
#include "test_common/has_ranges.hpp"
#include "test_common/network_result.hpp"
#include "test_integration/any_connection_fixture.hpp"

#ifdef BOOST_MYSQL_HAS_RANGES
#include <ranges>
#endif

using namespace boost::mysql;
using namespace boost::mysql::test;
using boost::test_tools::per_element;
namespace asio = boost::asio;

// Cover all possible execution requests for execute() and async_execute()

namespace {

BOOST_AUTO_TEST_SUITE(test_execution_requests)

BOOST_FIXTURE_TEST_CASE(query, any_connection_fixture)
{
    // Setup
    connect();
    results r;
    execution_state st;

    // execute
    conn.async_execute("SELECT 1", r, as_netresult).validate_no_error();
    BOOST_TEST(r.rows() == makerows(1, 1), per_element());

    // types convertible to string_view work
    conn.async_execute(std::string("SELECT 1"), r, as_netresult).validate_no_error();
    BOOST_TEST(r.rows() == makerows(1, 1), per_element());

    // spotcheck: start execution with a text query works
    conn.async_start_execution("SELECT 1", st, as_netresult).validate_no_error();
    auto rws = conn.async_read_some_rows(st, as_netresult).get();
    BOOST_TEST(rws == makerows(1, 1), per_element());
}

BOOST_FIXTURE_TEST_CASE(stmt_tuple, any_connection_fixture)
{
    // Also verifies that tuples correctly apply the writable field transformation
    // Setup
    connect();
    results r;
    execution_state st;
    auto stmt = conn.async_prepare_statement("SELECT ?, ?", as_netresult).get();
    BOOST_TEST_REQUIRE(stmt.num_params() == 2u);

    // execute
    conn.async_execute(stmt.bind("42", boost::optional<int>(13)), r, as_netresult).validate_no_error();
    BOOST_TEST(r.rows() == makerows(2, "42", 13), per_element());

    // references work
    std::string s = "abcdef";
    auto op = conn.async_execute(stmt.bind(std::ref(s), 21), r, asio::deferred);
    s = "opqrs";
    std::move(op)(as_netresult).validate_no_error();
    BOOST_TEST(r.rows() == makerows(2, "opqrs", 21), per_element());

    // spotcheck: start execution with tuples works
    conn.async_start_execution(stmt.bind("abc", boost::optional<int>()), st, as_netresult)
        .validate_no_error();
    auto rws = conn.async_read_some_rows(st, as_netresult).get();
    BOOST_TEST(rws == makerows(2, "abc", nullptr), per_element());

    // spotcheck: errors correctly detected
    conn.async_execute(stmt.bind("42"), r, as_netresult).validate_error(client_errc::wrong_num_params);

    // spotcheck: lvalues work
    const auto bound_stmt = stmt.bind("42", boost::optional<int>());
    conn.async_execute(bound_stmt, r, as_netresult).validate_no_error();
    BOOST_TEST(r.rows() == makerows(2, "42", nullptr), per_element());
}

BOOST_FIXTURE_TEST_CASE(stmt_range, any_connection_fixture)
{
    // Setup
    connect();
    results r;
    execution_state st;
    std::forward_list<field_view> params{field_view(42), field_view("abc")};
    auto stmt = conn.async_prepare_statement("SELECT ?, ?", as_netresult).get();
    BOOST_TEST_REQUIRE(stmt.num_params() == 2u);

    // execute
    conn.async_execute(stmt.bind(params.begin(), params.end()), r, as_netresult).validate_no_error();
    BOOST_TEST(r.rows() == makerows(2, 42, "abc"), per_element());

    // spotcheck: statement with ranges work with start execution
    conn.async_start_execution(stmt.bind(params.begin(), params.end()), st, as_netresult).validate_no_error();
    auto rws = conn.async_read_some_rows(st, as_netresult).get();
    BOOST_TEST(rws == makerows(2, 42, "abc"), per_element());

    // Regression check: executing with a type convertible (but not equal) to field_view works
    const std::vector<field> owning_params{field_view(50), field_view("luv")};
    conn.async_execute(stmt.bind(owning_params.begin(), owning_params.end()), r, as_netresult)
        .validate_no_error();
    BOOST_TEST(r.rows() == makerows(2, 50, "luv"), per_element());

    // Spotcheck: errors detected
    auto too_few_params = make_fv_arr(1);
    conn.async_execute(stmt.bind(too_few_params.begin(), too_few_params.end()), r, as_netresult)
        .validate_error(client_errc::wrong_num_params);

    // Spotchecks: lvalues work
    auto bound_stmt = stmt.bind(params.begin(), params.end());
    conn.async_execute(bound_stmt, r, as_netresult).validate_no_error();
    BOOST_TEST(r.rows() == makerows(2, 42, "abc"), per_element());
}

BOOST_FIXTURE_TEST_CASE(with_params_, any_connection_fixture)
{
    // Setup
    connect();
    results r;
    execution_state st;

    // execute
    conn.async_execute(with_params("SELECT {}, {}", 42, "abc"), r, as_netresult).validate_no_error();
    BOOST_TEST(r.rows() == makerows(2, 42, "abc"), per_element());

    // spotcheck: can be used with start_execution
    conn.async_start_execution(with_params("SELECT {}, {}", 42, "abc"), st, as_netresult).validate_no_error();
    auto rws = conn.async_read_some_rows(st, as_netresult).get();
    BOOST_TEST(rws == makerows(2, 42, "abc"), per_element());

    // references work
    std::string s = "abcdef";
    auto op = conn.async_execute(with_params("SELECT {}, {}", 42, std::ref(s)), r, asio::deferred);
    s = "opqrs";
    std::move(op)(as_netresult).validate_no_error();
    BOOST_TEST(r.rows() == makerows(2, 42, "opqrs"), per_element());

    // Queries without parameters work
    conn.async_execute(with_params("SELECT '{{}}'"), r, as_netresult).validate_no_error();
    BOOST_TEST(r.rows() == makerows(1, "{}"), per_element());

    // lvalues work
    const auto req = with_params("SELECT {}, {}", "42", boost::optional<int>(100));
    conn.async_execute(req, r, as_netresult).validate_no_error();
    BOOST_TEST(r.rows() == makerows(2, "42", 100), per_element());

#ifdef BOOST_MYSQL_HAS_RANGES
    // Regression check: mutable ranges work
    std::vector<int> nums{2, 5, 10, 20};
    auto is_even = [](int i) { return i % 2 == 0; };
    conn.async_execute(with_params("SELECT {}", std::ranges::views::filter(nums, is_even)), r, as_netresult)
        .validate_no_error();
    BOOST_TEST(r.rows() == makerows(3, 2, 10, 20), per_element());
#endif

    // Error: unformattable values
    conn.async_execute(with_params("SELECT {}", "bad\xffutf8"), r, as_netresult)
        .validate_error(client_errc::invalid_encoding);

    // Error: unknown charset
    conn.async_reset_connection(as_netresult).validate_no_error();
    conn.async_execute(with_params("SELECT {}", 42), r, as_netresult)
        .validate_error(client_errc::unknown_character_set);
}

// Spotcheck: with_params() is owning, and can be safely used
// together with sequence() in deferred ops without incurring in lifetime problems
BOOST_FIXTURE_TEST_CASE(with_params_sequence_deferred, any_connection_fixture)
{
    // Setup
    connect();
    results r;
    auto fn = [](int value, format_context_base& fmt_ctx) { fmt_ctx.append_value(value); };

    // Create a deferred op
    auto op = conn.async_execute(
        with_params("SELECT {}, {}", sequence(std::vector<int>{3, 4, 7}, fn, " + "), std::string(10, 'a')),
        r,
        asio::deferred
    );

    // Run it
    std::move(op)(as_netresult).validate_no_error();
    BOOST_TEST(r.rows() == makerows(2, 14, "aaaaaaaaaa"), per_element());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace