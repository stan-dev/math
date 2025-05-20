//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/results.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/io_context_fixture.hpp"
#include "test_common/network_result.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_ok_frame.hpp"
#include "test_unit/test_any_connection.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

BOOST_AUTO_TEST_SUITE(test_connection_use_after_move)

BOOST_FIXTURE_TEST_CASE(use_move_constructed_connection, io_context_fixture)
{
    // Construct connection
    auto conn = create_test_any_connection(ctx);

    // Use it
    get_stream(conn).add_bytes(create_ok_frame(1, ok_builder().build()));
    results result;
    conn.async_execute("SELECT * FROM myt", result, as_netresult).validate_no_error();

    // Move construct another connection from conn
    any_connection conn2(std::move(conn));

    // Using it works (no dangling pointer)
    get_stream(conn2).add_bytes(create_ok_frame(1, ok_builder().affected_rows(42).build()));
    conn2.async_execute("DELETE FROM myt", result, as_netresult).validate_no_error();
    BOOST_TEST(result.affected_rows() == 42u);
}

BOOST_FIXTURE_TEST_CASE(use_move_assigned_connection, io_context_fixture)
{
    // Construct two connections
    auto conn1 = create_test_any_connection(ctx);
    auto conn2 = create_test_any_connection(ctx);

    // Use them
    get_stream(conn1).add_bytes(create_ok_frame(1, ok_builder().build()));
    get_stream(conn2).add_bytes(create_ok_frame(1, ok_builder().build()));
    results result;
    conn1.async_execute("SELECT * FROM myt", result, as_netresult).validate_no_error();
    conn2.async_execute("SELECT * FROM myt", result, as_netresult).validate_no_error();

    // Move assign
    conn2 = std::move(conn1);

    // Using it works (no dangling pointer)
    get_stream(conn2).add_bytes(create_ok_frame(1, ok_builder().affected_rows(42).build()));
    conn2.async_execute("DELETE FROM myt", result, as_netresult).validate_no_error();
    BOOST_TEST(result.affected_rows() == 42u);
}

BOOST_AUTO_TEST_SUITE_END()
