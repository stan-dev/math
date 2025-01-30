//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/detail/config.hpp>

#ifdef BOOST_MYSQL_CXX14

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/static_execution_state.hpp>

#include <boost/core/span.hpp>
#include <boost/test/unit_test.hpp>

#include "test_common/io_context_fixture.hpp"
#include "test_common/network_result.hpp"
#include "test_unit/create_execution_processor.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_row_message.hpp"
#include "test_unit/test_any_connection.hpp"
#include "test_unit/test_stream.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;
using boost::span;

BOOST_AUTO_TEST_SUITE(test_read_some_rows_static)

using row1 = std::tuple<int, float>;
using row2 = std::tuple<double>;

using state_t = static_execution_state<row1, row1, row2, row1, row2>;

struct fixture : io_context_fixture
{
    state_t st;
    any_connection conn{create_test_any_connection(ctx)};
    std::array<row1, 3> storage1;
    std::array<row2, 3> storage2;

    span<row1> storage1_span() { return storage1; }
    span<row2> storage2_span() { return storage2; }

    test_stream& stream() noexcept { return get_stream(conn); }

    void add_ok() { ::add_ok(get_iface(st), ok_builder().more_results(true).build()); }

    void add_meta_row1()
    {
        add_meta(
            get_iface(st),
            {
                meta_builder().type(column_type::int_).nullable(false).build_coldef(),
                meta_builder().type(column_type::float_).nullable(false).build_coldef(),
            }
        );
    }

    void add_meta_row2()
    {
        add_meta(
            get_iface(st),
            {
                meta_builder().type(column_type::double_).nullable(false).build_coldef(),
            }
        );
    }
};

BOOST_FIXTURE_TEST_CASE(repeated_row_types, fixture)
{
    add_meta_row1();

    // 1st resultset: row1
    stream().add_bytes(create_text_row_message(0, 10, 4.2f)).add_bytes(create_text_row_message(1, 11, 4.3f));

    std::size_t num_rows = conn.async_read_some_rows(st, storage1_span(), as_netresult).get();
    BOOST_TEST_REQUIRE(num_rows == 2u);
    BOOST_TEST((storage1[0] == row1{10, 4.2f}));
    BOOST_TEST((storage1[1] == row1{11, 4.3f}));

    // Advance resultset
    add_ok();
    add_meta_row1();
    BOOST_TEST_REQUIRE(st.should_read_rows());

    // 2nd resultset: row1 again
    stream().add_bytes(create_text_row_message(2, 13, 0.2f));
    num_rows = conn.async_read_some_rows(st, storage1_span(), as_netresult).get();
    BOOST_TEST_REQUIRE(num_rows == 1u);
    BOOST_TEST((storage1[0] == row1{13, 0.2f}));

    // Advance resultset
    add_ok();
    add_meta_row2();
    BOOST_TEST_REQUIRE(st.should_read_rows());

    // 3rd resultset: row2
    stream().add_bytes(create_text_row_message(3, 9.1));
    num_rows = conn.async_read_some_rows(st, storage2_span(), as_netresult).get();
    BOOST_TEST_REQUIRE(num_rows == 1u);
    BOOST_TEST((storage2[0] == row2{9.1}));

    // Advance resultset
    add_ok();
    add_meta_row1();
    BOOST_TEST_REQUIRE(st.should_read_rows());

    // 4th resultset: row1
    stream().add_bytes(create_text_row_message(4, 43, 0.7f));
    num_rows = conn.async_read_some_rows(st, storage1_span(), as_netresult).get();
    BOOST_TEST_REQUIRE(num_rows == 1u);
    BOOST_TEST((storage1[0] == row1{43, 0.7f}));

    // Advance resultset
    add_ok();
    add_meta_row2();
    BOOST_TEST_REQUIRE(st.should_read_rows());

    // 5th resultset: row2
    stream().add_bytes(create_text_row_message(5, 99.9));
    num_rows = conn.async_read_some_rows(st, storage2_span(), as_netresult).get();
    BOOST_TEST_REQUIRE(num_rows == 1u);
    BOOST_TEST((storage2[0] == row2{99.9}));
}

BOOST_FIXTURE_TEST_CASE(error_row_type_mismatch, fixture)
{
    add_meta_row1();

    // 1st resultset: row1. Note that this will consume the message
    stream().add_bytes(create_text_row_message(0, 10, 4.2f));
    conn.async_read_some_rows(st, storage2_span(), as_netresult)
        .validate_error(client_errc::row_type_mismatch);

    // Advance resultset
    add_ok();
    add_meta_row1();
    add_ok();
    add_meta_row2();
    BOOST_TEST_REQUIRE(st.should_read_rows());

    // 3rd resultset: row2
    stream().add_bytes(create_text_row_message(1, 9.1));
    conn.async_read_some_rows(st, storage1_span(), as_netresult)
        .validate_error(client_errc::row_type_mismatch);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
