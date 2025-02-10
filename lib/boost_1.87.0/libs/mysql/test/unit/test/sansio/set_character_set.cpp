//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/character_set.hpp>
#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/impl/internal/sansio/set_character_set.hpp>

#include <boost/core/span.hpp>
#include <boost/test/unit_test.hpp>

#include "test_common/create_diagnostics.hpp"
#include "test_common/printing.hpp"
#include "test_unit/algo_test.hpp"
#include "test_unit/create_err.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_ok_frame.hpp"
#include "test_unit/create_query_frame.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;

BOOST_AUTO_TEST_SUITE(test_set_character_set)

//
// compose_set_names
//
BOOST_AUTO_TEST_CASE(compose_set_names_success)
{
    BOOST_TEST(detail::compose_set_names(utf8mb4_charset).value() == "SET NAMES 'utf8mb4'");
    BOOST_TEST(detail::compose_set_names(ascii_charset).value() == "SET NAMES 'ascii'");
}

BOOST_AUTO_TEST_CASE(compose_set_names_needs_escaping)
{
    // We don't create vulnerabilities when creating SET NAMES statements
    character_set mock_charset{"ab'cd\"e", utf8mb4_charset.next_char};
    BOOST_TEST(detail::compose_set_names(mock_charset).value() == "SET NAMES 'ab\\'cd\\\"e'");
}

BOOST_AUTO_TEST_CASE(compose_set_names_error)
{
    struct
    {
        string_view name;
        const char* charset_name;
    } test_cases[] = {
        {"utf8",     "test-\xc3\xb1"}, // anything non-ascii is rejected
        {"non_utf8", "test-\xb1-abc"},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.charset_name)
        {
            character_set charset{tc.charset_name, utf8mb4_charset.next_char};
            BOOST_TEST(detail::compose_set_names(charset).error() == client_errc::invalid_encoding);
        }
    }
}

//
// read_set_character_set_response_algo
//
struct read_response_fixture : algo_fixture_base
{
    detail::read_set_character_set_response_algo algo;

    // Clearing diagnostics is not this algorithm's responsibility
    read_response_fixture(character_set charset = utf8mb4_charset)
        : algo_fixture_base(diagnostics()), algo(diag, charset, 29)
    {
    }
};

BOOST_AUTO_TEST_CASE(read_response_success)
{
    // Setup
    read_response_fixture fix;

    // Run the algo
    algo_test().expect_read(create_ok_frame(29, ok_builder().build())).check(fix);

    // The charset was updated
    BOOST_TEST(fix.st.current_charset == utf8mb4_charset);
}

BOOST_AUTO_TEST_CASE(read_response_success_previous_charset)
{
    // Setup
    read_response_fixture fix;
    fix.st.current_charset = ascii_charset;

    // Run the algo
    algo_test().expect_read(create_ok_frame(29, ok_builder().build())).check(fix);

    // The charset was updated
    BOOST_TEST(fix.st.current_charset == utf8mb4_charset);
}

BOOST_AUTO_TEST_CASE(read_response_error_network)
{
    // This covers errors in read and write
    algo_test()
        .expect_read(create_ok_frame(29, ok_builder().build()))
        .check_network_errors<read_response_fixture>();
}

BOOST_AUTO_TEST_CASE(read_response_error_packet)
{
    // Setup
    read_response_fixture fix;

    // Run the algo
    algo_test()
        .expect_read(err_builder()
                         .seqnum(29)
                         .code(common_server_errc::er_unknown_character_set)
                         .message("Unknown charset")
                         .build_frame())
        .check(fix, common_server_errc::er_unknown_character_set, create_server_diag("Unknown charset"));

    // The current character set was not updated
    BOOST_TEST(fix.st.current_charset == character_set());
}

//
// set_character_set_algo
//
struct set_charset_fixture : algo_fixture_base
{
    detail::set_character_set_algo algo;

    set_charset_fixture(character_set charset = utf8mb4_charset) : algo(diag, {charset}) {}
    set_charset_fixture(std::size_t max_bufsize)
        : algo_fixture_base(max_bufsize), algo(diag, {utf8mb4_charset})
    {
    }
};

BOOST_AUTO_TEST_CASE(set_charset_success)
{
    // Setup
    set_charset_fixture fix;

    // Run the algo
    algo_test()
        .expect_write(create_query_frame(0, "SET NAMES 'utf8mb4'"))
        .expect_read(create_ok_frame(1, ok_builder().build()))
        .check(fix);

    // The charset was updated
    BOOST_TEST(fix.st.current_charset == utf8mb4_charset);
}

// Ensure we don't create vulnerabilities when composing SET NAMES
BOOST_AUTO_TEST_CASE(set_charset_name_needs_escaping)
{
    // Setup
    set_charset_fixture fix({"lat'in\\", utf8mb4_charset.next_char});

    // Run the algo
    algo_test()
        .expect_write(create_query_frame(0, "SET NAMES 'lat\\'in\\\\'"))
        .expect_read(create_ok_frame(1, ok_builder().build()))
        .check(fix);
}

BOOST_AUTO_TEST_CASE(set_charset_error_composing_request)
{
    // Setup
    set_charset_fixture fix({"lat\xc3\xadn", utf8mb4_charset.next_char});
    fix.st.current_charset = ascii_charset;

    // Run the algo
    algo_test().check(fix, client_errc::invalid_encoding);

    // The current character set was not updated
    BOOST_TEST(fix.st.current_charset == ascii_charset);
}

BOOST_AUTO_TEST_CASE(set_charset_error_network)
{
    // This covers errors in reading the response
    algo_test()
        .expect_write(create_query_frame(0, "SET NAMES 'utf8mb4'"))
        .expect_read(create_ok_frame(1, ok_builder().build()))
        .check_network_errors<set_charset_fixture>();
}

BOOST_AUTO_TEST_CASE(set_charset_error_max_buffer_size)
{
    // Setup
    set_charset_fixture fix(16u);

    // Run the algo
    algo_test().check(fix, client_errc::max_buffer_size_exceeded);
}

BOOST_AUTO_TEST_SUITE_END()
