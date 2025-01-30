//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/blob.hpp>
#include <boost/mysql/column_type.hpp>
#include <boost/mysql/connection.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/mysql_collations.hpp>
#include <boost/mysql/pipeline.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/rows_view.hpp>
#include <boost/mysql/statement.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/detail/access.hpp>
#include <boost/mysql/detail/engine.hpp>
#include <boost/mysql/detail/engine_impl.hpp>
#include <boost/mysql/detail/engine_stream_adaptor.hpp>

#include <boost/asio/deferred.hpp>
#include <boost/asio/error.hpp>
#include <boost/test/unit_test.hpp>

#include <memory>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/buffer_concat.hpp"
#include "test_common/io_context_fixture.hpp"
#include "test_common/netfun_maker.hpp"
#include "test_common/network_result.hpp"
#include "test_common/printing.hpp"
#include "test_unit/create_coldef_frame.hpp"
#include "test_unit/create_execution_processor.hpp"
#include "test_unit/create_frame.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_ok_frame.hpp"
#include "test_unit/create_row_message.hpp"
#include "test_unit/create_statement.hpp"
#include "test_unit/fail_count.hpp"
#include "test_unit/test_any_connection.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;
namespace asio = boost::asio;

BOOST_AUTO_TEST_SUITE(test_misc)

// spotcheck for the dynamic interface
// Verifies that execute (dynamic interface) works when rows come in separate batches
// This is testing the interaction between the network algorithm and results
BOOST_FIXTURE_TEST_CASE(execute_multiple_batches, io_context_fixture)
{
    // Setup
    auto conn = create_test_any_connection(ctx);
    results result;

    // Message sequence (each on its own read)
    get_stream(conn)
        .add_bytes(create_frame(1, {0x02}))  // OK, 2 columns
        .add_break()
        .add_bytes(create_coldef_frame(2, meta_builder().type(column_type::varchar).build_coldef()))  // meta
        .add_break()
        .add_bytes(create_coldef_frame(
            3,
            meta_builder().type(column_type::blob).collation_id(mysql_collations::binary).build_coldef()
        ))  // meta
        .add_break()
        .add_bytes(create_text_row_message(4, "abcd", makebv("\0\1\0")))  // row 1
        .add_break()
        .add_bytes(create_text_row_message(5, "defghi", makebv("\3\4\3\0")))  // row 2
        .add_break()
        .add_bytes(create_eof_frame(6, ok_builder().affected_rows(10u).info("1st").more_results(true).build())
        )
        .add_break()
        .add_bytes(create_ok_frame(7, ok_builder().affected_rows(20u).info("2nd").more_results(true).build()))
        .add_break()
        .add_bytes(create_frame(8, {0x01}))  // OK, 1 metadata
        .add_break()
        .add_bytes(create_coldef_frame(9, meta_builder().type(column_type::varchar).build_coldef()))  // meta
        .add_break()
        .add_bytes(create_text_row_message(10, "ab"))  // row 1
        .add_break()
        .add_bytes(create_eof_frame(11, ok_builder().affected_rows(30u).info("3rd").build()));

    // Call the function
    conn.execute("abc", result);

    // We've written the query request
    auto expected_msg = create_frame(0, {0x03, 0x61, 0x62, 0x63});  // ASCII "abc" (plus length)
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(get_stream(conn).bytes_written(), expected_msg);

    // We've populated the results
    BOOST_TEST_REQUIRE(result.size() == 3u);
    BOOST_TEST(result[0].affected_rows() == 10u);
    BOOST_TEST(result[0].info() == "1st");
    BOOST_TEST(result[0].rows() == makerows(2, "abcd", makebv("\0\1\0"), "defghi", makebv("\3\4\3\0")));
    BOOST_TEST(result[1].affected_rows() == 20u);
    BOOST_TEST(result[1].info() == "2nd");
    BOOST_TEST(result[1].rows() == rows());
    BOOST_TEST(result[2].affected_rows() == 30u);
    BOOST_TEST(result[2].info() == "3rd");
    BOOST_TEST(result[2].rows() == makerows(1, "ab"));
}

// Regression check: async_close_statement doesn't require the passed-in statement to be alive
// when used with deferred tokens.
BOOST_FIXTURE_TEST_CASE(async_close_statement_handle_deferred_tokens, io_context_fixture)
{
    // Setup
    auto conn = create_test_any_connection(ctx);
    std::unique_ptr<statement> stmt{new statement(statement_builder().id(3).build())};
    get_stream(conn).add_bytes(create_ok_frame(1, ok_builder().build()));

    // Deferred op
    auto op = conn.async_close_statement(*stmt, asio::deferred);

    // Invalidate the original variable
    stmt.reset();

    // Run the operation
    std::move(op)(as_netresult).validate_no_error();

    // verify that the op had the intended effects
    const auto expected_message = concat(
        create_frame(0, {0x19, 0x03, 0x00, 0x00, 0x00}),
        create_frame(0, {0x0e})
    );
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(get_stream(conn).bytes_written(), expected_message);
}

// Regression check: when there is a network error, sync functions
// returning a value fail with an assertion
BOOST_AUTO_TEST_CASE(net_error_prepare_statement)
{
    using netmaker_stmt = netfun_maker<statement, any_connection, string_view>;
    struct
    {
        const char* name;
        netmaker_stmt::signature prepare_statement;
    } fns[] = {
        {"sync",  netmaker_stmt::sync_errc(&any_connection::prepare_statement)       },
        {"async", netmaker_stmt::async_diag(&any_connection::async_prepare_statement)},
    };

    for (const auto& fn : fns)
    {
        BOOST_TEST_CONTEXT(fn.name)
        {
            // Setup
            io_context_fixture fix;
            auto conn = create_test_any_connection(fix.ctx);
            get_stream(conn).set_fail_count(fail_count(0, boost::asio::error::connection_reset));

            fn.prepare_statement(conn, "SELECT 1").validate_error(boost::asio::error::connection_reset);
        }
    }
}

BOOST_AUTO_TEST_CASE(net_error_read_some_rows)
{
    using netmaker_stmt = netfun_maker<rows_view, any_connection, execution_state&>;
    struct
    {
        const char* name;
        netmaker_stmt::signature read_some_rows;
    } fns[] = {
        {"sync",  netmaker_stmt::sync_errc(&any_connection::read_some_rows)       },
        {"async", netmaker_stmt::async_diag(&any_connection::async_read_some_rows)},
    };

    for (const auto& fn : fns)
    {
        BOOST_TEST_CONTEXT(fn.name)
        {
            // Setup
            io_context_fixture fix;
            auto conn = create_test_any_connection(fix.ctx);
            get_stream(conn).set_fail_count(fail_count(0, boost::asio::error::connection_reset));
            execution_state st;
            add_meta(get_iface(st), {column_type::bigint});

            fn.read_some_rows(conn, st).validate_error(boost::asio::error::connection_reset);
        }
    }
}

BOOST_AUTO_TEST_CASE(net_error_void_signature)
{
    using netmaker_execute = netfun_maker<void, any_connection, const string_view&, results&>;
    struct
    {
        const char* name;
        netmaker_execute::signature execute;
    } fns[] = {
        {"sync",  netmaker_execute::sync_errc(&any_connection::execute)       },
        {"async", netmaker_execute::async_diag(&any_connection::async_execute)},
    };

    for (const auto& fn : fns)
    {
        BOOST_TEST_CONTEXT(fn.name)
        {
            // Setup
            io_context_fixture fix;
            auto conn = create_test_any_connection(fix.ctx);
            get_stream(conn).set_fail_count(fail_count(0, boost::asio::error::connection_reset));
            results r;

            fn.execute(conn, "SELECT 1", r).validate_error(boost::asio::error::connection_reset);
        }
    }
}

// empty pipelines complete immediately, posting adequately
BOOST_FIXTURE_TEST_CASE(empty_pipeline, io_context_fixture)
{
    // Setup
    auto conn = create_test_any_connection(ctx);
    pipeline_request req;
    std::vector<stage_response> res;

    // Run it. It should complete immediately, posting to the correct executor (verified by the testing
    // infrastructure)
    conn.async_run_pipeline(req, res, as_netresult).validate_no_error();
    BOOST_TEST(res.size() == 0u);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(get_stream(conn).bytes_written(), blob{});
}

// fatal errors in pipelines behave correctly
BOOST_FIXTURE_TEST_CASE(pipeline_fatal_error, io_context_fixture)
{
    // Setup
    auto conn = create_test_any_connection(ctx);
    pipeline_request req;
    std::vector<stage_response> res;
    req.add_execute("SELECT 1").add_execute("SELECT 2");

    // The first read will fail
    get_stream(conn).set_fail_count(fail_count(1, boost::asio::error::network_reset));

    // Run it
    conn.async_run_pipeline(req, res, as_netresult).validate_error(boost::asio::error::network_reset);

    // Validate the results
    BOOST_TEST(res.size() == 2u);
    BOOST_TEST(res[0].error() == boost::asio::error::network_reset);
    BOOST_TEST(res[0].diag() == diagnostics());
    BOOST_TEST(res[1].error() == boost::asio::error::network_reset);
    BOOST_TEST(res[1].diag() == diagnostics());
}

BOOST_AUTO_TEST_SUITE_END()
