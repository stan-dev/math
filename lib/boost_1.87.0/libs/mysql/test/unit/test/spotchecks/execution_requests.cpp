//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

// Since integration tests can't reliably test multifunction operations
// that span over multiple messages, we test the complete multifn fllow in this unit tests.

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/character_set.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/rows_view.hpp>
#include <boost/mysql/string_view.hpp>
#include <boost/mysql/tcp.hpp>
#include <boost/mysql/with_params.hpp>

#include <boost/mysql/detail/access.hpp>

#include <boost/asio/deferred.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/buffer_concat.hpp"
#include "test_common/io_context_fixture.hpp"
#include "test_common/network_result.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_ok_frame.hpp"
#include "test_unit/create_query_frame.hpp"
#include "test_unit/create_statement.hpp"
#include "test_unit/test_any_connection.hpp"
#include "test_unit/test_stream.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;
namespace asio = boost::asio;

BOOST_AUTO_TEST_SUITE(test_execution_requests)

// The execution request is forwarded correctly, and non-const objects work correctly.
// This is relevant for with_params involving ranges::views::filter, for instance
// Note: netmakers intentionally not used here
BOOST_AUTO_TEST_SUITE(forwarding_constness)

struct fixture : io_context_fixture
{
    any_connection conn{create_test_any_connection(ctx)};
    results result;
    execution_state st;

    fixture() { get_stream(conn).add_bytes(create_ok_frame(1, ok_builder().build())); }

    void validate_bytes_written()
    {
        BOOST_MYSQL_ASSERT_BUFFER_EQUALS(
            get_stream(conn).bytes_written(),
            create_query_frame(0, "SELECT 'abc'")
        );
    }

    struct request
    {
        string_view value;
        operator string_view() { return value; }  // Intentionally non-const
    };
    static request make_request() { return {"SELECT 'abc'"}; }
};

BOOST_FIXTURE_TEST_CASE(execute_sync_errc, fixture)
{
    // Setup
    error_code ec;
    diagnostics diag;

    // Run
    conn.execute(make_request(), result, ec, diag);

    // Validate
    BOOST_TEST(ec == error_code());
    BOOST_TEST(diag == diagnostics());
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(execute_sync_exc, fixture)
{
    conn.execute(make_request(), result);
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(execute_async_diag, fixture)
{
    diagnostics diag;
    conn.async_execute(make_request(), result, diag, as_netresult).validate_no_error();
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(execute_async_nodiag, fixture)
{
    conn.async_execute(make_request(), result, as_netresult).validate_no_error();
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(execute_async_deferred, fixture)
{
    // Spotcheck: deferred works correctly
    auto op = conn.async_execute(make_request(), result, asio::deferred);
    std::move(op)(as_netresult).validate_no_error();
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(start_execution_sync_errc, fixture)
{
    // Setup
    error_code ec;
    diagnostics diag;

    // Run
    conn.start_execution(make_request(), st, ec, diag);

    // Validate
    BOOST_TEST(ec == error_code());
    BOOST_TEST(diag == diagnostics());
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(start_execution_sync_exc, fixture)
{
    conn.start_execution(make_request(), st);
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(start_execution_async_diag, fixture)
{
    diagnostics diag;
    conn.async_start_execution(make_request(), st, diag, as_netresult).validate_no_error();
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(start_execution_async_nodiag, fixture)
{
    conn.async_start_execution(make_request(), st, as_netresult).validate_no_error();
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(start_execution_async_deferred, fixture)
{
    // Spotcheck: deferred works correctly
    auto op = conn.async_start_execution(make_request(), st, asio::deferred);
    std::move(op)(as_netresult).validate_no_error();
    validate_bytes_written();
}

// Spotcheck: old connection doesn't generate compile errors.
// Intentionally not run
void old_connection()
{
    asio::io_context ctx;
    tcp_connection conn(ctx);
    auto req = fixture::make_request();
    results result;
    execution_state st;
    error_code ec;
    diagnostics diag;

    conn.execute(req, result, ec, diag);
    conn.execute(req, result);
    conn.async_execute(req, result, asio::detached);
    conn.async_execute(req, result, diag, asio::detached);
    conn.async_execute(req, result, diag, asio::deferred)(asio::detached);
    conn.start_execution(req, st, ec, diag);
    conn.start_execution(req, st);
    conn.async_start_execution(req, st, asio::detached);
    conn.async_start_execution(req, st, diag, asio::detached);
    conn.async_start_execution(req, st, diag, asio::deferred)(asio::detached);
}

BOOST_AUTO_TEST_SUITE_END()

// The execution request is forwarded with the correct value category, and is moved if required.
// Move-only objects work correctly
BOOST_AUTO_TEST_SUITE(forwarding_rvalues)

struct fixture : io_context_fixture
{
    any_connection conn{create_test_any_connection(ctx)};
    results result;
    execution_state st;

    fixture() { get_stream(conn).add_bytes(create_ok_frame(1, ok_builder().build())); }

    void validate_bytes_written()
    {
        BOOST_MYSQL_ASSERT_BUFFER_EQUALS(
            get_stream(conn).bytes_written(),
            create_query_frame(0, "SELECT 'abcd'")
        );
    }

    // An execution request that doesn't support copies
    struct request
    {
        std::unique_ptr<char[]> buff;
        std::size_t size;

        request(string_view v) : buff(new char[v.size()]), size(v.size())
        {
            std::copy(v.begin(), v.end(), buff.get());
        }

        operator string_view() const { return {buff.get(), size}; }
    };
    static request make_request() { return {"SELECT 'abcd'"}; }
};

BOOST_FIXTURE_TEST_CASE(execute_sync_errc, fixture)
{
    // Setup
    error_code ec;
    diagnostics diag;

    // Execute
    conn.execute(make_request(), result, ec, diag);

    // Check
    BOOST_TEST(ec == error_code());
    BOOST_TEST(diag == diagnostics());
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(execute_sync_exc, fixture)
{
    conn.execute(make_request(), result);
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(execute_async_diag, fixture)
{
    diagnostics diag;
    conn.async_execute(make_request(), result, diag, as_netresult).validate_no_error();
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(execute_async_nodiag, fixture)
{
    conn.async_execute(make_request(), result, as_netresult).validate_no_error();
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(execute_async_deferred, fixture)
{
    auto op = conn.async_execute(make_request(), result, asio::deferred);
    std::move(op)(as_netresult).validate_no_error();
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(start_execution_sync_errc, fixture)
{
    // Setup
    error_code ec;
    diagnostics diag;

    // Execute
    conn.start_execution(make_request(), st, ec, diag);

    // Check
    BOOST_TEST(ec == error_code());
    BOOST_TEST(diag == diagnostics());
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(start_execution_sync_exc, fixture)
{
    conn.start_execution(make_request(), st);
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(start_execution_async_diag, fixture)
{
    diagnostics diag;
    conn.async_start_execution(make_request(), st, diag, as_netresult).validate_no_error();
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(start_execution_async_nodiag, fixture)
{
    conn.async_start_execution(make_request(), st, as_netresult).validate_no_error();
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(start_execution_async_deferred, fixture)
{
    auto op = conn.async_start_execution(make_request(), st, asio::deferred);
    std::move(op)(as_netresult).validate_no_error();
    validate_bytes_written();
}
// Spotcheck: old connection doesn't generate compile errors.
// Intentionally not run
void old_connection()
{
    asio::io_context ctx;
    tcp_connection conn(ctx);
    results result;
    execution_state st;
    error_code ec;
    diagnostics diag;

    conn.execute(fixture::make_request(), result, ec, diag);
    conn.execute(fixture::make_request(), result);
    conn.async_execute(fixture::make_request(), result, asio::detached);
    conn.async_execute(fixture::make_request(), result, diag, asio::detached);
    conn.async_execute(fixture::make_request(), result, diag, asio::deferred)(asio::detached);
    conn.start_execution(fixture::make_request(), st, ec, diag);
    conn.start_execution(fixture::make_request(), st);
    conn.async_start_execution(fixture::make_request(), st, asio::detached);
    conn.async_start_execution(fixture::make_request(), st, diag, asio::detached);
    conn.async_start_execution(fixture::make_request(), st, diag, asio::deferred)(asio::detached);
}

BOOST_AUTO_TEST_SUITE_END()

// The execution request is forwarded correctly, with the correct value category.
// Lvalues are not moved
BOOST_AUTO_TEST_SUITE(forwarding_lvalues)

struct fixture : io_context_fixture
{
    // A request where we can detect moved-from state.
    // std::string move-constructor doesn't offer guarantees about moved-from objects,
    // while std::vector does
    struct vector_request
    {
        std::vector<char> buff;

        vector_request(string_view v) : buff(v.begin(), v.end()) {}

        operator string_view() const { return {buff.data(), buff.size()}; }
    };

    any_connection conn{create_test_any_connection(ctx)};
    results result;
    execution_state st;
    vector_request req{"SELECT 'abcd'"};

    fixture() { get_stream(conn).add_bytes(create_ok_frame(1, ok_builder().build())); }

    void validate_bytes_written()
    {
        BOOST_MYSQL_ASSERT_BUFFER_EQUALS(
            get_stream(conn).bytes_written(),
            create_query_frame(0, "SELECT 'abcd'")
        );
        BOOST_TEST(static_cast<string_view>(req) == "SELECT 'abcd'");
    }
};

BOOST_FIXTURE_TEST_CASE(execute_sync_errc, fixture)
{
    // Setup
    error_code ec;
    diagnostics diag;

    // Execute
    conn.execute(req, result, ec, diag);

    // Validate
    BOOST_TEST(ec == error_code());
    BOOST_TEST(diag == diagnostics());
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(execute_sync_exc, fixture)
{
    conn.execute(req, result);
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(execute_async_diag, fixture)
{
    diagnostics diag;
    conn.async_execute(req, result, diag, as_netresult).validate_no_error();
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(execute_async_nodiag, fixture)
{
    conn.async_execute(req, result, as_netresult).validate_no_error();
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(execute_async_deferred, fixture)
{
    auto op = conn.async_execute(req, result, asio::deferred);
    std::move(op)(as_netresult).validate_no_error();
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(start_execution_sync_errc, fixture)
{
    // Setup
    error_code ec;
    diagnostics diag;

    // Execute
    conn.start_execution(req, st, ec, diag);

    // Check
    BOOST_TEST(ec == error_code());
    BOOST_TEST(diag == diagnostics());
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(start_execution_sync_exc, fixture)
{
    conn.start_execution(req, st);
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(start_execution_async_diag, fixture)
{
    diagnostics diag;
    conn.async_start_execution(req, st, diag, as_netresult).validate_no_error();
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(start_execution_async_nodiag, fixture)
{
    conn.async_start_execution(req, st, as_netresult).validate_no_error();
    validate_bytes_written();
}

BOOST_FIXTURE_TEST_CASE(start_execution_async_deferred, fixture)
{
    auto op = conn.async_start_execution(req, st, asio::deferred);
    std::move(op)(as_netresult).validate_no_error();
    validate_bytes_written();
}

// Spotcheck: old connection doesn't generate compile errors.
// Intentionally not run
void old_connection()
{
    asio::io_context ctx;
    tcp_connection conn(ctx);
    results result;
    execution_state st;
    error_code ec;
    diagnostics diag;
    const fixture::vector_request req{"SELECT 'abc'"};

    conn.execute(req, result, ec, diag);
    conn.execute(req, result);
    conn.async_execute(req, result, asio::detached);
    conn.async_execute(req, result, diag, asio::detached);
    conn.async_execute(req, result, diag, asio::deferred)(asio::detached);
    conn.start_execution(req, st, ec, diag);
    conn.start_execution(req, st);
    conn.async_start_execution(req, st, asio::detached);
    conn.async_start_execution(req, st, diag, asio::detached);
    conn.async_start_execution(req, st, diag, asio::deferred)(asio::detached);
}

BOOST_AUTO_TEST_SUITE_END()

// Deferred tokens appropriately decay-copy lvalues
BOOST_AUTO_TEST_SUITE(deferred_tokens_lvalues)

BOOST_FIXTURE_TEST_CASE(execute, io_context_fixture)
{
    // Setup
    auto conn = create_test_any_connection(ctx);
    results result;
    get_stream(conn).add_bytes(create_ok_frame(1, ok_builder().build()));
    std::string req(128, 'a');

    // Create a deferred op
    auto op = conn.async_execute(req, result, asio::deferred);

    // Mutate the argument
    std::fill(req.begin(), req.end(), 'b');

    // Initiate
    std::move(op)(as_netresult).validate_no_error();

    // We wrote the initial value
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(
        get_stream(conn).bytes_written(),
        create_query_frame(0, std::string(128, 'a'))
    );
}

BOOST_FIXTURE_TEST_CASE(start_execution, io_context_fixture)
{
    // Setup
    auto conn = create_test_any_connection(ctx);
    execution_state st;
    get_stream(conn).add_bytes(create_ok_frame(1, ok_builder().build()));
    std::string req(128, 'a');

    // Create a deferred op
    auto op = conn.async_start_execution(req, st, asio::deferred);

    // Mutate the argument
    std::fill(req.begin(), req.end(), 'b');

    // Initiate
    std::move(op)(as_netresult).validate_no_error();

    // We wrote the initial value
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(
        get_stream(conn).bytes_written(),
        create_query_frame(0, std::string(128, 'a'))
    );
}

BOOST_AUTO_TEST_SUITE_END()

// Spotcheck: the types returned by with_params are correct
BOOST_AUTO_TEST_CASE(with_params_types)
{
    {
        // const references
        const std::string s = "abc";
        auto p1 = with_params("SELECT {}", s);
        static_assert(std::is_same<decltype(p1), with_params_t<std::string>>::value, "");
    }
    {
        // references
        std::string s = "abc";
        auto p = with_params("SELECT {}", s);
        static_assert(std::is_same<decltype(p), with_params_t<std::string>>::value, "");
    }
    {
        // rvalue references
        std::string s = "abc";
        auto p = with_params("SELECT {}", std::move(s));
        static_assert(std::is_same<decltype(p), with_params_t<std::string>>::value, "");
    }
    {
        // pure rvalues
        auto p = with_params("SELECT {}", std::string());
        static_assert(std::is_same<decltype(p), with_params_t<std::string>>::value, "");
    }
    {
        // std::ref
        std::string s = "abc";
        auto p = with_params("SELECT {}", std::ref(s));
        static_assert(std::is_same<decltype(p), with_params_t<std::string&>>::value, "");
    }
    {
        // std::ref (const)
        const std::string s = "abc";
        auto p = with_params("SELECT {}", std::ref(s));
        static_assert(std::is_same<decltype(p), with_params_t<const std::string&>>::value, "");
    }
}

// Regression test: async_execute() doesn't cause side effects in the initiation
BOOST_FIXTURE_TEST_CASE(async_execute_side_effects_in_initiation, io_context_fixture)
{
    auto conn = create_test_any_connection(ctx);
    results result1, result2;

    // Resultsets will be complete as soon as a message is read
    get_stream(conn)
        .add_bytes(create_ok_frame(1, ok_builder().affected_rows(2).build()))
        .add_bytes(create_ok_frame(1, ok_builder().affected_rows(1).build()));

    // Create two queries as deferred objects, but don't run them yet
    auto q1 = conn.async_execute("Q1", result1, asio::deferred);
    auto q2 = conn.async_execute("Q2", result2, asio::deferred);

    // Run them in reverse order
    std::move(q2)(as_netresult).validate_no_error();
    std::move(q1)(as_netresult).validate_no_error();

    // Check that we wrote Q2's message first, then Q1's
    auto expected = concat(
        create_frame(0, {0x03, 0x51, 0x32}),  // query request Q2
        create_frame(0, {0x03, 0x51, 0x31})   // query request Q1
    );
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(get_stream(conn).bytes_written(), expected);

    // Check that the results got the right ok_packets
    BOOST_TEST(result2.affected_rows() == 2u);
    BOOST_TEST(result1.affected_rows() == 1u);
}

// Regression test: bound statements correctly store statement handle and params
// when used with deferred tokens
BOOST_FIXTURE_TEST_CASE(async_execute_deferred_lifetimes, io_context_fixture)
{
    // Setup
    constexpr std::uint8_t expected_msg[] = {
        0x15, 0x00, 0x00, 0x00, 0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00,
        0x00, 0x02, 0x01, 0xfe, 0x00, 0x06, 0x00, 0x04, 0x74, 0x65, 0x73, 0x74,
    };
    results result;
    auto conn = create_test_any_connection(ctx);
    get_stream(conn).add_bytes(create_ok_frame(1, ok_builder().build()));

    // Create a bound statement on the heap. This helps tooling detect memory errors
    using bound_stmt_t = bound_statement_tuple<std::tuple<std::string, std::nullptr_t>>;
    std::unique_ptr<bound_stmt_t> stmt_ptr{
        new bound_stmt_t{statement_builder().id(1).num_params(2).build().bind(std::string("test"), nullptr)}
    };

    // Deferred op
    auto op = conn.async_execute(*stmt_ptr, result, asio::deferred);

    // Free the statement
    stmt_ptr.reset();

    // Actually run the op
    std::move(op)(as_netresult).validate_no_error();

    // verify that the op had the intended effects
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(get_stream(conn).bytes_written(), expected_msg);
}

BOOST_AUTO_TEST_SUITE_END()
