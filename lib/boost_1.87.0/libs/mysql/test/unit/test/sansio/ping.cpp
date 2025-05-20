//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>

#include <boost/mysql/impl/internal/sansio/connection_state_data.hpp>
#include <boost/mysql/impl/internal/sansio/ping.hpp>
#include <boost/mysql/impl/internal/sansio/run_pipeline.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/printing.hpp"
#include "test_unit/algo_test.hpp"
#include "test_unit/create_err.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_ok_frame.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;

BOOST_AUTO_TEST_SUITE(test_ping)

//
// read_ping_response_algo
//
struct read_response_fixture : algo_fixture_base
{
    detail::read_ping_response_algo algo{diag, 57};

    // Clearing diagnostics is not this algorithm's responsibility
    read_response_fixture() : algo_fixture_base(diagnostics()) {}
};

BOOST_AUTO_TEST_CASE(read_response_success)
{
    // Setup
    read_response_fixture fix;

    // Run the test
    algo_test()
        .expect_read(create_ok_frame(57, ok_builder().build()))  // OK response
        .check(fix);

    // The OK packet was processed correctly
    BOOST_TEST(fix.st.backslash_escapes);
}

BOOST_AUTO_TEST_CASE(read_response_success_no_backslash_escapes)
{
    // Setup
    read_response_fixture fix;

    // Run the test
    algo_test()
        .expect_read(create_ok_frame(57, ok_builder().no_backslash_escapes(true).build()))  // OK response
        .check(fix);

    // The OK packet was processed correctly
    BOOST_TEST(!fix.st.backslash_escapes);
}

BOOST_AUTO_TEST_CASE(read_response_error_network)
{
    algo_test()
        .expect_read(create_ok_frame(57, ok_builder().build()))
        .check_network_errors<read_response_fixture>();
}

BOOST_AUTO_TEST_CASE(read_response_error_packet)
{
    // Setup
    read_response_fixture fix;

    // Run the test
    algo_test()
        .expect_read(err_builder()
                         .seqnum(57)
                         .code(common_server_errc::er_bad_db_error)
                         .message("my_message")
                         .build_frame())  // Error response
        .check(fix, common_server_errc::er_bad_db_error, create_server_diag("my_message"));
}

//
// setup_ping_pipeline: running a pipeline with these parameters
// has the intended effect
//

struct ping_fixture : algo_fixture_base
{
    detail::run_pipeline_algo algo{diag, detail::setup_ping_pipeline(st)};
};

BOOST_AUTO_TEST_CASE(ping_success)
{
    // Setup
    ping_fixture fix;

    // Run the test
    algo_test()
        .expect_write({0x01, 0x00, 0x00, 0x00, 0x0e})           // ping request
        .expect_read(create_ok_frame(1, ok_builder().build()))  // OK response
        .check(fix);

    // The OK packet was processed correctly
    BOOST_TEST(fix.st.backslash_escapes);
}

BOOST_AUTO_TEST_CASE(ping_error_network)
{
    // Check for net errors for each read/write
    algo_test()
        .expect_write({0x01, 0x00, 0x00, 0x00, 0x0e})
        .expect_read(create_ok_frame(1, ok_builder().build()))
        .check_network_errors<ping_fixture>();
}

BOOST_AUTO_TEST_CASE(ping_error_response)
{
    // Setup
    ping_fixture fix;

    // Run the test
    algo_test()
        .expect_write({0x01, 0x00, 0x00, 0x00, 0x0e})  // Ping request
        .expect_read(err_builder()
                         .seqnum(1)
                         .code(common_server_errc::er_bad_db_error)
                         .message("my_message")
                         .build_frame())  // Error response
        .check(fix, common_server_errc::er_bad_db_error, create_server_diag("my_message"));
}

BOOST_AUTO_TEST_SUITE_END()
