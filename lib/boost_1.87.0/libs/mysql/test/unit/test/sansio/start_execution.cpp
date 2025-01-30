//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/character_set.hpp>
#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/column_type.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/metadata_mode.hpp>

#include <boost/mysql/detail/any_execution_request.hpp>
#include <boost/mysql/detail/resultset_encoding.hpp>

#include <boost/mysql/impl/internal/sansio/connection_state_data.hpp>
#include <boost/mysql/impl/internal/sansio/start_execution.hpp>

#include <boost/test/unit_test.hpp>

#include <array>
#include <string>

#include "test_common/check_meta.hpp"
#include "test_common/create_basic.hpp"
#include "test_unit/algo_test.hpp"
#include "test_unit/create_coldef_frame.hpp"
#include "test_unit/create_frame.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_ok_frame.hpp"
#include "test_unit/create_query_frame.hpp"
#include "test_unit/mock_execution_processor.hpp"
#include "test_unit/printing.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;
using boost::mysql::detail::any_execution_request;
using boost::mysql::detail::resultset_encoding;

BOOST_AUTO_TEST_SUITE(test_start_execution)

struct fixture : algo_fixture_base
{
    mock_execution_processor proc;
    detail::start_execution_algo algo;

    fixture(
        any_execution_request req = any_execution_request("SELECT 1"),
        std::size_t max_bufsize = default_max_buffsize
    )
        : algo_fixture_base(max_bufsize), algo(diag, {req, &proc})
    {
    }
};

BOOST_AUTO_TEST_CASE(text_query)
{
    // Setup
    fixture fix(any_execution_request("SELECT 1"));

    // Run the algo
    algo_test()
        .expect_write(create_frame(0, {0x03, 0x53, 0x45, 0x4c, 0x45, 0x43, 0x54, 0x20, 0x31}))
        .expect_read(create_frame(1, {0x01}))
        .expect_read(create_coldef_frame(2, meta_builder().type(column_type::varchar).build_coldef()))
        .check(fix);

    // Verify
    BOOST_TEST(fix.proc.encoding() == resultset_encoding::text);
    BOOST_TEST(fix.proc.sequence_number() == 3u);
    BOOST_TEST(fix.proc.is_reading_rows());
    check_meta(fix.proc.meta(), {column_type::varchar});
    fix.proc.num_calls().reset(1).on_num_meta(1).on_meta(1).validate();
}

BOOST_AUTO_TEST_CASE(stmt_success)
{
    // Setup
    const auto params = make_fv_arr("test", nullptr);
    fixture fix(any_execution_request({std::uint32_t(1u), std::uint16_t(2u), params}));

    // Run the algo
    algo_test()
        .expect_write(create_frame(
            0,
            {
                0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x02,
                0x01, 0xfe, 0x00, 0x06, 0x00, 0x04, 0x74, 0x65, 0x73, 0x74,
            }
        ))
        .expect_read(create_frame(1, {0x01}))
        .expect_read(create_coldef_frame(2, meta_builder().type(column_type::varchar).build_coldef()))
        .check(fix);

    // Verify
    BOOST_TEST(fix.proc.encoding() == resultset_encoding::binary);
    BOOST_TEST(fix.proc.sequence_number() == 3u);
    BOOST_TEST(fix.proc.is_reading_rows());
    check_meta(fix.proc.meta(), {column_type::varchar});
    fix.proc.num_calls().reset(1).on_num_meta(1).on_meta(1).validate();
}

BOOST_AUTO_TEST_CASE(stmt_error_num_params)
{
    // Setup
    const auto params = make_fv_arr("test", nullptr, 42);  // too many params
    fixture fix(any_execution_request({std::uint32_t(1u), std::uint16_t(2u), params}));

    // Run the algo. Nothing should be written to the server
    algo_test().check(fix, client_errc::wrong_num_params);
}

BOOST_AUTO_TEST_CASE(with_params_success)
{
    // Setup
    const std::array<format_arg, 2> args{
        {{"", "abc"}, {"", 42}}
    };
    fixture fix(any_execution_request({"SELECT {}, {}", args}));
    fix.st.current_charset = utf8mb4_charset;

    // Run the algo
    algo_test()
        .expect_write(create_query_frame(0, "SELECT 'abc', 42"))
        .expect_read(create_ok_frame(1, ok_builder().build()))
        .check(fix);

    // Verify
    BOOST_TEST(fix.proc.encoding() == resultset_encoding::text);
    BOOST_TEST(fix.proc.sequence_number() == 2u);
    BOOST_TEST(fix.proc.is_complete());
    fix.proc.num_calls().reset(1).on_head_ok_packet(1).validate();
}

BOOST_AUTO_TEST_CASE(with_params_error_unknown_charset)
{
    // Setup
    const std::array<format_arg, 2> args{
        {{"", "abc"}, {"", 42}}
    };
    fixture fix(any_execution_request({"SELECT {}, {}", args}));
    fix.st.current_charset = {};

    // The algo fails immediately
    algo_test().check(fix, client_errc::unknown_character_set);
}

BOOST_AUTO_TEST_CASE(with_params_error_formatting)
{
    // Setup
    const std::array<format_arg, 1> args{{{"", "abc"}}};
    fixture fix(any_execution_request({"SELECT {}, {}", args}));
    fix.st.current_charset = utf8mb4_charset;

    // The algo fails immediately
    algo_test().check(fix, client_errc::format_arg_not_found);
}

// This covers errors in both writing the request and calling read_resultset_head
BOOST_AUTO_TEST_CASE(error_network_error)
{
    // This will test for errors writing the execution request
    // and reading the response and metadata (thus, calling read_resultset_head)
    algo_test()
        .expect_write(create_frame(0, {0x03, 0x53, 0x45, 0x4c, 0x45, 0x43, 0x54, 0x20, 0x31}))
        .expect_read(create_frame(1, {0x01}))
        .expect_read(create_coldef_frame(2, meta_builder().type(column_type::varchar).build_coldef()))
        .check_network_errors<fixture>();
}

BOOST_AUTO_TEST_CASE(error_max_buffer_size)
{
    // Setup
    std::string query(512, 'a');
    fixture fix(any_execution_request(query), 512u);

    // Run the algo
    algo_test().check(fix, client_errc::max_buffer_size_exceeded);
}

BOOST_AUTO_TEST_SUITE_END()
