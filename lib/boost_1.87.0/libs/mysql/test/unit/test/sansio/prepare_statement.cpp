//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/statement.hpp>

#include <boost/mysql/impl/internal/protocol/impl/protocol_types.hpp>
#include <boost/mysql/impl/internal/protocol/impl/serialization_context.hpp>
#include <boost/mysql/impl/internal/sansio/prepare_statement.hpp>

#include <boost/test/unit_test.hpp>

#include "test_unit/algo_test.hpp"
#include "test_unit/create_coldef_frame.hpp"
#include "test_unit/create_err.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_prepare_statement_response.hpp"
#include "test_unit/create_query_frame.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;

BOOST_AUTO_TEST_SUITE(test_prepare_statement)

//
// read_prepare_statement_response_algo
//
struct read_response_fixture : algo_fixture_base
{
    detail::read_prepare_statement_response_algo algo{diag, 19};

    // Clearing diagnostics is not this algorithm's responsibility
    read_response_fixture() : algo_fixture_base(diagnostics()) {}

    statement result() const { return algo.result(st); }
};

BOOST_AUTO_TEST_CASE(read_response_success)
{
    // Setup
    read_response_fixture fix;

    // Run the algo
    algo_test()
        .expect_read(prepare_stmt_response_builder().seqnum(19).id(1).num_columns(1).num_params(2).build())
        .expect_read(create_coldef_frame(20, meta_builder().name("abc").build_coldef()))
        .expect_read(create_coldef_frame(21, meta_builder().name("other").build_coldef()))
        .expect_read(create_coldef_frame(22, meta_builder().name("final").build_coldef()))
        .check(fix);

    // The statement was created successfully
    auto stmt = fix.result();
    BOOST_TEST(stmt.id() == 1u);
    BOOST_TEST(stmt.num_params() == 2u);
}

BOOST_AUTO_TEST_CASE(read_response_success_0cols)
{
    // Setup
    read_response_fixture fix;

    // Run the algo
    algo_test()
        .expect_read(prepare_stmt_response_builder().seqnum(19).id(5).num_columns(0).num_params(1).build())
        .expect_read(create_coldef_frame(20, meta_builder().name("abc").build_coldef()))
        .check(fix);

    // The statement was created successfully
    auto stmt = fix.result();
    BOOST_TEST(stmt.id() == 5u);
    BOOST_TEST(stmt.num_params() == 1u);
}

BOOST_AUTO_TEST_CASE(read_response_success_0params)
{
    // Setup
    read_response_fixture fix;

    // Run the algo
    algo_test()
        .expect_read(prepare_stmt_response_builder().seqnum(19).id(214).num_columns(2).num_params(0).build())
        .expect_read(create_coldef_frame(20, meta_builder().name("abc").build_coldef()))
        .expect_read(create_coldef_frame(21, meta_builder().name("defff").build_coldef()))
        .check(fix);

    // The statement was created successfully
    auto stmt = fix.result();
    BOOST_TEST(stmt.id() == 214u);
    BOOST_TEST(stmt.num_params() == 0u);
}

BOOST_AUTO_TEST_CASE(read_response_success_0cols_0params)
{
    // Setup
    read_response_fixture fix;

    // Run the algo
    algo_test()
        .expect_read(prepare_stmt_response_builder().seqnum(19).id(98).num_columns(0).num_params(0).build())
        .check(fix);

    // The statement was created successfully
    auto stmt = fix.result();
    BOOST_TEST(stmt.id() == 98u);
    BOOST_TEST(stmt.num_params() == 0u);
}

// A statement_id == 0 doesn't cause trouble
BOOST_AUTO_TEST_CASE(read_response_success_0id)
{
    // Setup
    read_response_fixture fix;

    // Run the algo
    algo_test()
        .expect_read(prepare_stmt_response_builder().seqnum(19).id(0).num_columns(0).num_params(1).build())
        .expect_read(create_coldef_frame(20, meta_builder().name("abc").build_coldef()))
        .check(fix);

    // The statement was created successfully
    auto stmt = fix.result();
    BOOST_TEST(stmt.id() == 0u);
    BOOST_TEST(stmt.num_params() == 1u);
}

BOOST_AUTO_TEST_CASE(read_response_error_network)
{
    algo_test()
        .expect_read(prepare_stmt_response_builder().seqnum(19).id(1).num_columns(1).num_params(2).build())
        .expect_read(create_coldef_frame(20, meta_builder().name("abc").build_coldef()))
        .expect_read(create_coldef_frame(21, meta_builder().name("other").build_coldef()))
        .expect_read(create_coldef_frame(22, meta_builder().name("final").build_coldef()))
        .check_network_errors<read_response_fixture>();
}

BOOST_AUTO_TEST_CASE(read_response_error_packet)
{
    // Setup
    read_response_fixture fix;

    // Run the algo
    algo_test()
        .expect_read(err_builder()
                         .seqnum(19)
                         .code(common_server_errc::er_bad_db_error)
                         .message("my_message")
                         .build_frame())
        .check(fix, common_server_errc::er_bad_db_error, create_server_diag("my_message"));
}

//
// prepare_statement_algo
//
struct prepare_fixture : algo_fixture_base
{
    detail::prepare_statement_algo algo{diag, {"SELECT 1"}};

    prepare_fixture() = default;
    prepare_fixture(std::size_t max_bufsize) : algo_fixture_base(max_bufsize) {}

    statement result() const { return algo.result(st); }
};

BOOST_AUTO_TEST_CASE(prepare_success)
{
    // Setup
    prepare_fixture fix;

    // Run the algo
    algo_test()
        .expect_write(create_prepare_statement_frame(0, "SELECT 1"))
        .expect_read(prepare_stmt_response_builder().seqnum(1).id(29).num_columns(0).num_params(2).build())
        .expect_read(create_coldef_frame(2, meta_builder().name("abc").build_coldef()))
        .expect_read(create_coldef_frame(3, meta_builder().name("other").build_coldef()))
        .check(fix);

    // The statement was created successfully
    auto stmt = fix.result();
    BOOST_TEST(stmt.id() == 29u);
    BOOST_TEST(stmt.num_params() == 2u);
}

// Spotcheck: an error while reading the response is propagated correctly
BOOST_AUTO_TEST_CASE(prepare_error_packet)
{
    // Setup
    prepare_fixture fix;

    // Run the algo
    algo_test()
        .expect_write(create_prepare_statement_frame(0, "SELECT 1"))
        .expect_read(err_builder()
                         .seqnum(1)
                         .code(common_server_errc::er_bad_db_error)
                         .message("my_message")
                         .build_frame())
        .check(fix, common_server_errc::er_bad_db_error, create_server_diag("my_message"));
}

BOOST_AUTO_TEST_CASE(prepare_network_error)
{
    // This covers errors in the request and the response
    algo_test()
        .expect_write(create_prepare_statement_frame(0, "SELECT 1"))
        .expect_read(prepare_stmt_response_builder().seqnum(1).id(29).num_columns(0).num_params(1).build())
        .expect_read(create_coldef_frame(2, meta_builder().name("abc").build_coldef()))
        .check_network_errors<prepare_fixture>();
}

BOOST_AUTO_TEST_CASE(prepare_error_max_buffer_size)
{
    // Setup
    prepare_fixture fix(10u);

    // Run the algo
    algo_test().check(fix, client_errc::max_buffer_size_exceeded);
}

BOOST_AUTO_TEST_SUITE_END()
