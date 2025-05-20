//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/detail/config.hpp>

#ifdef BOOST_MYSQL_CXX14

#include <boost/mysql/column_type.hpp>
#include <boost/mysql/static_execution_state.hpp>

#include <boost/describe/class.hpp>
#include <boost/describe/operators.hpp>
#include <boost/test/unit_test.hpp>

#include "test_common/check_meta.hpp"
#include "test_unit/create_execution_processor.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

BOOST_AUTO_TEST_SUITE(test_static_execution_state)

struct row1
{
    std::string f1;
    float f2;
};
BOOST_DESCRIBE_STRUCT(row1, (), (f1, f2))

using boost::describe::operators::operator==;
using boost::describe::operators::operator<<;

using row2 = std::tuple<double>;
using empty = std::tuple<>;

using execst_t = static_execution_state<row1, row2, empty>;

// The functionality has been tested in static_execution_state_impl already.
// Just spotchecks here
BOOST_AUTO_TEST_CASE(spotchecks)
{
    execst_t st;
    auto& impl = get_iface(st);
    diagnostics diag;

    // Initial
    BOOST_TEST(st.should_start_op());
    BOOST_TEST(!st.should_read_head());
    BOOST_TEST(!st.should_read_rows());
    BOOST_TEST(!st.complete());
    BOOST_TEST(st.meta().size() == 0u);

    // Meta
    add_meta(
        impl,
        {
            meta_builder().type(column_type::varchar).nullable(false).name("f1").build_coldef(),
            meta_builder().type(column_type::float_).nullable(false).name("f2").build_coldef(),
        }
    );
    BOOST_TEST(!st.should_start_op());
    BOOST_TEST(!st.should_read_head());
    BOOST_TEST(st.should_read_rows());
    BOOST_TEST(!st.complete());
    check_meta(st.meta(), {column_type::varchar, column_type::float_});

    // First resultset OK
    add_ok(
        impl,
        ok_builder()
            .affected_rows(1)
            .last_insert_id(2)
            .warnings(4)
            .info("abc")
            .more_results(true)
            .out_params(true)
            .build()
    );
    BOOST_TEST(!st.should_start_op());
    BOOST_TEST(st.should_read_head());
    BOOST_TEST(!st.should_read_rows());
    BOOST_TEST(!st.complete());
    check_meta(st.meta(), {column_type::varchar, column_type::float_});
    BOOST_TEST(st.affected_rows() == 1u);
    BOOST_TEST(st.last_insert_id() == 2u);
    BOOST_TEST(st.warning_count() == 4u);
    BOOST_TEST(st.info() == "abc");
    BOOST_TEST(st.is_out_params());

    // Second resultset meta
    add_meta(impl, {meta_builder().type(column_type::double_).nullable(false).build_coldef()});
    BOOST_TEST(!st.should_start_op());
    BOOST_TEST(!st.should_read_head());
    BOOST_TEST(st.should_read_rows());
    BOOST_TEST(!st.complete());
    check_meta(st.meta(), {column_type::double_});

    // Second resultset OK
    add_ok(
        impl,
        ok_builder().affected_rows(4).last_insert_id(5).warnings(6).info("2nd").more_results(true).build()
    );
    BOOST_TEST(!st.should_start_op());
    BOOST_TEST(st.should_read_head());
    BOOST_TEST(!st.should_read_rows());
    BOOST_TEST(!st.complete());
    check_meta(st.meta(), {column_type::double_});
    BOOST_TEST(st.affected_rows() == 4u);
    BOOST_TEST(st.last_insert_id() == 5u);
    BOOST_TEST(st.warning_count() == 6u);
    BOOST_TEST(st.info() == "2nd");
    BOOST_TEST(!st.is_out_params());

    // Third, empty resultset
    add_ok(impl, ok_builder().info("3rd").build());
    BOOST_TEST(!st.should_start_op());
    BOOST_TEST(!st.should_read_head());
    BOOST_TEST(!st.should_read_rows());
    BOOST_TEST(st.complete());
    BOOST_TEST(st.meta().empty());
    BOOST_TEST(st.affected_rows() == 0u);
    BOOST_TEST(st.last_insert_id() == 0u);
    BOOST_TEST(st.warning_count() == 0u);
    BOOST_TEST(st.info() == "3rd");
    BOOST_TEST(!st.is_out_params());
}

std::unique_ptr<execst_t> create_heap_state()
{
    std::unique_ptr<execst_t> res{new execst_t};
    add_meta(
        get_iface(*res),
        {
            meta_builder().type(column_type::varchar).nullable(false).name("f1").build_coldef(),
            meta_builder().type(column_type::float_).nullable(false).name("f2").build_coldef(),
        }
    );
    add_ok(
        get_iface(*res),
        ok_builder().affected_rows(1).last_insert_id(2).warnings(4).info("1st").more_results(true).build()
    );
    return res;
}

// Verify that the lifetime guarantees we make are correct
BOOST_AUTO_TEST_CASE(move_constructor)
{
    // Having this in heap helps detect lifetime issues
    auto st = create_heap_state();

    // Obtain references
    auto meta = st->meta();
    auto info = st->info();

    // Move construct
    execst_t st2(std::move(*st));
    st.reset();

    // Make sure that views are still valid
    check_meta(meta, {column_type::varchar, column_type::float_});
    BOOST_TEST(info == "1st");

    // The new object holds the same data
    BOOST_TEST_REQUIRE(st2.should_read_head());
    check_meta(st2.meta(), {column_type::varchar, column_type::float_});
    BOOST_TEST(st2.info() == "1st");
}

BOOST_AUTO_TEST_CASE(move_assignment)
{
    // Having this in heap helps detect lifetime issues
    auto st = create_heap_state();

    // Obtain references
    auto meta = st->meta();
    auto info = st->info();

    // Move assign
    execst_t st2;
    st2 = std::move(*st);
    st.reset();

    // Make sure that views are still valid
    check_meta(meta, {column_type::varchar, column_type::float_});
    BOOST_TEST(info == "1st");

    // The new object holds the same data
    BOOST_TEST_REQUIRE(st2.should_read_head());
    check_meta(st2.meta(), {column_type::varchar, column_type::float_});
    BOOST_TEST(st2.info() == "1st");
}

BOOST_AUTO_TEST_SUITE_END()

#endif
