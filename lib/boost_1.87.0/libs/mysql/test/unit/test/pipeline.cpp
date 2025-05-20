//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/character_set.hpp>
#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/pipeline.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/detail/access.hpp>
#include <boost/mysql/detail/pipeline.hpp>
#include <boost/mysql/detail/resultset_encoding.hpp>

#include <boost/core/ignore_unused.hpp>
#include <boost/core/span.hpp>
#include <boost/optional/optional.hpp>
#include <boost/test/unit_test.hpp>

#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/buffer_concat.hpp"
#include "test_common/create_basic.hpp"
#include "test_common/create_diagnostics.hpp"
#include "test_common/printing.hpp"
#include "test_unit/create_execution_processor.hpp"
#include "test_unit/create_frame.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_query_frame.hpp"
#include "test_unit/create_statement.hpp"
#include "test_unit/printing.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;
using boost::test_tools::per_element;
using detail::pipeline_request_stage;
using detail::pipeline_stage_kind;
using detail::resultset_encoding;

BOOST_AUTO_TEST_SUITE(test_pipeline)

// Validates the exception thrown when adding a statement execution with wrong number of params
static auto stmt_exc_validator = [](const std::invalid_argument& exc) {
    BOOST_TEST(
        string_view(exc.what()) == "Wrong number of actual parameters supplied to a prepared statement"
    );
    return true;
};

BOOST_AUTO_TEST_SUITE(stage_response_)

BOOST_AUTO_TEST_CASE(default_ctor)
{
    // Construct
    stage_response r;

    // Contains an empty error
    BOOST_TEST(!r.has_results());
    BOOST_TEST(!r.has_statement());
    BOOST_TEST(r.error() == error_code());
    BOOST_TEST(r.diag() == diagnostics());
    BOOST_TEST(std::move(r).diag() == diagnostics());
}

BOOST_AUTO_TEST_CASE(underlying_error)
{
    // Setup
    stage_response r;
    detail::access::get_impl(r).set_error(client_errc::invalid_encoding, create_server_diag("my_message"));

    // Check
    BOOST_TEST(!r.has_results());
    BOOST_TEST(!r.has_statement());
    BOOST_TEST(r.error() == client_errc::invalid_encoding);
    BOOST_TEST(r.diag() == create_server_diag("my_message"));
    BOOST_TEST(std::move(r).diag() == create_server_diag("my_message"));
}

BOOST_AUTO_TEST_CASE(underlying_statement)
{
    // Setup
    stage_response r;
    detail::access::get_impl(r).set_result(statement_builder().id(3).build());

    // Check
    BOOST_TEST(!r.has_results());
    BOOST_TEST(r.has_statement());
    BOOST_TEST(r.as_statement().id() == 3u);
    BOOST_TEST(r.get_statement().id() == 3u);

    // error(), diag() can be called and return empty objects
    BOOST_TEST(r.error() == error_code());
    BOOST_TEST(r.diag() == diagnostics());
    BOOST_TEST(std::move(r).diag() == diagnostics());
}

BOOST_AUTO_TEST_CASE(underlying_results)
{
    // Setup
    stage_response r;
    detail::access::get_impl(r).emplace_results();
    add_ok(detail::access::get_impl(r).get_processor(), ok_builder().info("some_info").build());

    // Check
    BOOST_TEST(r.has_results());
    BOOST_TEST(!r.has_statement());
    BOOST_TEST(r.get_results().info() == "some_info");
    BOOST_TEST(r.as_results().info() == "some_info");

    // Rvalue reference accessors work
    results&& ref1 = std::move(r).get_results();
    results&& ref2 = std::move(r).as_results();
    boost::ignore_unused(ref1);
    boost::ignore_unused(ref2);

    // error(), diag() can be called and return empty objects
    BOOST_TEST(r.error() == error_code());
    BOOST_TEST(r.diag() == diagnostics());
    BOOST_TEST(std::move(r).diag() == diagnostics());
}

BOOST_AUTO_TEST_CASE(as_results_error)
{
    // Empty error
    stage_response r;
    BOOST_CHECK_THROW(r.as_results(), std::invalid_argument);
    BOOST_CHECK_THROW(std::move(r).as_results(), std::invalid_argument);

    // Non-empty error
    detail::access::get_impl(r).set_error(client_errc::extra_bytes, create_client_diag("my_msg"));
    BOOST_CHECK_THROW(r.as_results(), std::invalid_argument);
    BOOST_CHECK_THROW(std::move(r).as_results(), std::invalid_argument);

    // Statement
    detail::access::get_impl(r).set_result(statement_builder().build());
    BOOST_CHECK_THROW(r.as_results(), std::invalid_argument);
    BOOST_CHECK_THROW(std::move(r).as_results(), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(as_statement_error)
{
    // Empty error
    stage_response r;
    BOOST_CHECK_THROW(r.as_statement(), std::invalid_argument);

    // Non-empty error
    detail::access::get_impl(r).set_error(client_errc::extra_bytes, create_client_diag("my_msg"));
    BOOST_CHECK_THROW(r.as_statement(), std::invalid_argument);

    // results
    detail::access::get_impl(r).emplace_results();
    BOOST_CHECK_THROW(r.as_statement(), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(change_type)
{
    stage_response r;

    // Set results
    detail::access::get_impl(r).emplace_results();
    BOOST_TEST(r.has_results());

    // Set an error
    detail::access::get_impl(r).set_error(client_errc::extra_bytes, create_client_diag("abc"));
    BOOST_TEST(!r.has_results());
    BOOST_TEST(r.error() == client_errc::extra_bytes);
    BOOST_TEST(r.diag() == create_client_diag("abc"));

    // Reset the error
    detail::access::get_impl(r).emplace_error();
    BOOST_TEST(r.error() == error_code());
    BOOST_TEST(r.diag() == diagnostics());

    // Set a statement
    detail::access::get_impl(r).set_result(statement_builder().build());
    BOOST_TEST(r.has_statement());

    // Set results again
    detail::access::get_impl(r).emplace_results();
    BOOST_TEST(r.has_results());
    BOOST_TEST(!r.has_statement());
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(pipeline_request_)

// Helper to check the contents of a pipeline_request
void check_pipeline(
    const pipeline_request& req,
    const std::vector<std::uint8_t>& expected_buffer,
    boost::span<const detail::pipeline_request_stage> expected_stages
)
{
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(detail::access::get_impl(req).buffer_, expected_buffer);
    BOOST_TEST(detail::access::get_impl(req).stages_ == expected_stages, per_element());
}

// Helper for pipelines with a single request
void check_pipeline_single(
    const pipeline_request& req,
    const std::vector<std::uint8_t>& expected_buffer,
    detail::pipeline_request_stage expected_stage
)
{
    check_pipeline(req, expected_buffer, {&expected_stage, 1});
}

// Text query
BOOST_AUTO_TEST_CASE(add_execute_text_query)
{
    pipeline_request req;
    req.add_execute("SELECT 1");
    check_pipeline_single(
        req,
        create_query_frame(0, "SELECT 1"),
        {pipeline_stage_kind::execute, 1u, resultset_encoding::text}
    );
}

// Statement, passing parameters as individual arguments
BOOST_AUTO_TEST_CASE(add_execute_statement)
{
    pipeline_request req;
    req.add_execute(statement_builder().id(2).num_params(3).build(), 42, "abc", nullptr);
    check_pipeline_single(
        req,
        {0x1e, 0x00, 0x00, 0x00, 0x17, 0x02, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00,
         0x00, 0x00, 0x04, 0x01, 0x08, 0x00, 0xfe, 0x00, 0x06, 0x00, 0x2a, 0x00,
         0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0x61, 0x62, 0x63},
        {pipeline_stage_kind::execute, 1u, resultset_encoding::binary}
    );
}

BOOST_AUTO_TEST_CASE(add_execute_statement_writable_fields)
{
    // We run the required writable field transformations
    pipeline_request req;
    req.add_execute(
        statement_builder().id(2).num_params(3).build(),
        boost::optional<int>(42),
        "abc",
        boost::optional<int>()
    );
    check_pipeline_single(
        req,
        {0x1e, 0x00, 0x00, 0x00, 0x17, 0x02, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00,
         0x00, 0x00, 0x04, 0x01, 0x08, 0x00, 0xfe, 0x00, 0x06, 0x00, 0x2a, 0x00,
         0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0x61, 0x62, 0x63},
        {pipeline_stage_kind::execute, 1u, resultset_encoding::binary}
    );
}

BOOST_AUTO_TEST_CASE(execute_statement_no_params)
{
    pipeline_request req;
    req.add_execute(statement_builder().id(2).num_params(0).build());
    check_pipeline_single(
        req,
        {0x0a, 0x00, 0x00, 0x00, 0x17, 0x02, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00},
        {pipeline_stage_kind::execute, 1u, resultset_encoding::binary}
    );
}

BOOST_AUTO_TEST_CASE(add_execute_statement_too_few_params)
{
    pipeline_request req;
    BOOST_CHECK_EXCEPTION(
        req.add_execute(statement_builder().num_params(2).build(), 10),
        std::invalid_argument,
        stmt_exc_validator
    );
    check_pipeline(req, {}, {});  // Request unmodified
}

BOOST_AUTO_TEST_CASE(add_execute_statement_too_many_params)
{
    pipeline_request req;
    BOOST_CHECK_EXCEPTION(
        req.add_execute(statement_builder().num_params(2).build(), 10, 20, 30),
        std::invalid_argument,
        stmt_exc_validator
    );
    check_pipeline(req, {}, {});  // Request unmodified
}

// Statement, passing parameters as a range
BOOST_AUTO_TEST_CASE(add_execute_statement_range)
{
    pipeline_request req;
    req.add_execute_range(statement_builder().id(2).num_params(3).build(), make_fv_arr(42, "abc", nullptr));
    check_pipeline_single(
        req,
        {0x1e, 0x00, 0x00, 0x00, 0x17, 0x02, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00,
         0x00, 0x00, 0x04, 0x01, 0x08, 0x00, 0xfe, 0x00, 0x06, 0x00, 0x2a, 0x00,
         0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0x61, 0x62, 0x63},
        {pipeline_stage_kind::execute, 1u, resultset_encoding::binary}
    );
}

BOOST_AUTO_TEST_CASE(add_execute_statement_range_too_few_params)
{
    pipeline_request req;
    BOOST_CHECK_EXCEPTION(
        req.add_execute_range(statement_builder().num_params(2).build(), make_fv_arr(42)),
        std::invalid_argument,
        stmt_exc_validator
    );
    check_pipeline(req, {}, {});  // Request unmodified
}

BOOST_AUTO_TEST_CASE(add_execute_statement_range_too_many_params)
{
    pipeline_request req;
    BOOST_CHECK_EXCEPTION(
        req.add_execute_range(statement_builder().num_params(2).build(), make_fv_arr(42, nullptr, "abc")),
        std::invalid_argument,
        stmt_exc_validator
    );
    check_pipeline(req, {}, {});  // Request unmodified
}

// prepare statement
BOOST_AUTO_TEST_CASE(add_prepare_statement)
{
    pipeline_request req;
    req.add_prepare_statement("SELECT 1");
    check_pipeline_single(
        req,
        create_prepare_statement_frame(0, "SELECT 1"),
        {pipeline_stage_kind::prepare_statement, 1u, {}}
    );
}

BOOST_AUTO_TEST_CASE(add_prepare_statement_empty)
{
    // Empty prepare statements don't produce problems
    pipeline_request req;
    req.add_prepare_statement("");
    check_pipeline_single(
        req,
        create_prepare_statement_frame(0, ""),
        {pipeline_stage_kind::prepare_statement, 1u, {}}
    );
}

// close statement
BOOST_AUTO_TEST_CASE(add_close_statement)
{
    pipeline_request req;
    req.add_close_statement(statement_builder().id(3).num_params(1).build());
    check_pipeline_single(
        req,
        create_frame(0, {0x19, 0x03, 0x00, 0x00, 0x00}),
        {pipeline_stage_kind::close_statement, 1u, {}}
    );
}

// reset connection
BOOST_AUTO_TEST_CASE(add_reset_connection)
{
    pipeline_request req;
    req.add_reset_connection();
    check_pipeline_single(req, create_frame(0, {0x1f}), {pipeline_stage_kind::reset_connection, 1u, {}});
}

// set character set
BOOST_AUTO_TEST_CASE(add_set_character_set)
{
    pipeline_request req;
    req.add_set_character_set(utf8mb4_charset);
    check_pipeline_single(
        req,
        create_query_frame(0, "SET NAMES 'utf8mb4'"),
        {pipeline_stage_kind::set_character_set, 1u, utf8mb4_charset}
    );
}

BOOST_AUTO_TEST_CASE(add_set_character_set_escapes)
{
    // We don't create SQL injection vulnerabilities while composing SET NAMES
    character_set charset{"inj'ection", utf8mb4_charset.next_char};
    pipeline_request req;
    req.add_set_character_set(charset);
    check_pipeline_single(
        req,
        create_query_frame(0, "SET NAMES 'inj\\'ection'"),
        {pipeline_stage_kind::set_character_set, 1u, charset}
    );
}

BOOST_AUTO_TEST_CASE(add_set_character_set_error)
{
    // If a character set name that can't be securely escaped gets passed, we throw
    pipeline_request req;
    BOOST_CHECK_EXCEPTION(
        req.add_set_character_set(character_set{"bad\xff", utf8mb4_charset.next_char}),
        std::invalid_argument,
        [](const std::invalid_argument& exc) {
            BOOST_TEST(string_view(exc.what()) == "Invalid character set name");
            return true;
        }
    );
    check_pipeline(req, {}, {});  // request unmodified
}

// Several stages
BOOST_AUTO_TEST_CASE(add_incrementally)
{
    // Default ctor creates an empty request
    pipeline_request req;
    check_pipeline(req, {}, {});

    // Add a reset connection stage
    req.add_reset_connection();
    std::vector<pipeline_request_stage> expected_stages{
        {pipeline_stage_kind::reset_connection, 1u, {}}
    };
    check_pipeline(req, create_frame(0, {0x1f}), expected_stages);

    // Add an execution stage
    req.add_execute("SELECT 1");
    expected_stages = {
        {pipeline_stage_kind::reset_connection, 1u, {}                      },
        {pipeline_stage_kind::execute,          1u, resultset_encoding::text}
    };
    check_pipeline(req, concat(create_frame(0, {0x1f}), create_query_frame(0, "SELECT 1")), expected_stages);
}

BOOST_AUTO_TEST_CASE(all_stage_kinds)
{
    // Setup
    pipeline_request req;

    // Add stages
    req.add_reset_connection()
        .add_execute("SELECT 1")
        .add_prepare_statement("SELECT ?")
        .add_set_character_set(utf8mb4_charset)
        .add_close_statement(statement_builder().id(8).build());

    // Check
    auto expected_buffer = buffer_builder()
                               .add(create_frame(0, {0x1f}))
                               .add(create_query_frame(0, "SELECT 1"))
                               .add(create_prepare_statement_frame(0, "SELECT ?"))
                               .add(create_query_frame(0, "SET NAMES 'utf8mb4'"))
                               .add(create_frame(0, {0x19, 0x08, 0x00, 0x00, 0x00}))
                               .build();
    const std::array<pipeline_request_stage, 5> expected_stages{
        {
         {pipeline_stage_kind::reset_connection, 1u, {}},
         {pipeline_stage_kind::execute, 1u, resultset_encoding::text},
         {pipeline_stage_kind::prepare_statement, 1u, {}},
         {pipeline_stage_kind::set_character_set, 1u, utf8mb4_charset},
         {pipeline_stage_kind::close_statement, 1u, {}},
         }
    };
    check_pipeline(req, expected_buffer, expected_stages);
}

BOOST_AUTO_TEST_CASE(clear)
{
    // Create a pipeline request with some steps
    pipeline_request req;
    req.add_reset_connection().add_set_character_set(utf8mb4_charset).add_execute("SELECT 1");

    // Clear the pipeline
    req.clear();
    check_pipeline(req, {}, {});

    // Add some stages again
    req.add_execute("abc").add_close_statement(statement_builder().id(7).build());

    // Check
    const std::array<pipeline_request_stage, 2> expected_stages{
        {
         {pipeline_stage_kind::execute, 1u, resultset_encoding::text},
         {pipeline_stage_kind::close_statement, 1u, {}},
         }
    };
    check_pipeline(
        req,
        concat(create_query_frame(0, "abc"), create_frame(0, {0x19, 0x07, 0x00, 0x00, 0x00})),
        expected_stages
    );
}

BOOST_AUTO_TEST_CASE(clear_empty)
{
    // Clearing an empty pipeline is a no-op
    pipeline_request req;
    req.clear();
    check_pipeline(req, {}, {});
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
