//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/detail/config.hpp>

#ifdef BOOST_MYSQL_CXX14

#include <boost/mysql/column_type.hpp>
#include <boost/mysql/static_results.hpp>

#include <boost/describe/class.hpp>
#include <boost/describe/operators.hpp>

#include <tuple>

#include "test_common/check_meta.hpp"
#include "test_unit/create_execution_processor.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

BOOST_AUTO_TEST_SUITE(test_static_results)

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

using results_t = static_results<row1, row2, empty>;

results_t create_initial_results()
{
    results_t res;
    exec_access(get_iface(res))
        .meta({
            meta_builder().type(column_type::varchar).nullable(false).name("f1").build_coldef(),
            meta_builder().type(column_type::float_).nullable(false).name("f2").build_coldef(),
        })
        .row("abc", 10.1f)
        .row("def", 20.0f)
        .ok(ok_builder().affected_rows(1).last_insert_id(2).warnings(3).info("1st").more_results(true).build()
        )
        .meta({meta_builder().type(column_type::double_).nullable(false).build_coldef()})
        .row(42.0)
        .ok(ok_builder().affected_rows(4).last_insert_id(5).warnings(6).info("2nd").more_results(true).build()
        )
        .ok(ok_builder().info("3rd").build());
    return res;
}

std::unique_ptr<results_t> create_heap_results()
{
    return std::unique_ptr<results_t>(new results_t(create_initial_results()));
}

BOOST_AUTO_TEST_CASE(has_value)
{
    // Default construction
    results_t result;
    BOOST_TEST_REQUIRE(!result.has_value());

    // With value
    result = create_initial_results();
    BOOST_TEST_REQUIRE(result.has_value());
}

BOOST_AUTO_TEST_CASE(accessors)
{
    auto result = create_initial_results();

    // Rows
    BOOST_TEST_REQUIRE(result.rows<0>().size() == 2u);
    BOOST_TEST((result.rows<0>()[0] == row1{"abc", 10.1f}));
    BOOST_TEST((result.rows<0>()[1] == row1{"def", 20.0f}));
    BOOST_TEST_REQUIRE(result.rows<1>().size() == 1u);
    BOOST_TEST((result.rows<1>()[0] == row2{42.0}));
    BOOST_TEST((result.rows<2>().size() == 0u));

    // Meta
    check_meta(result.meta<0>(), {column_type::varchar, column_type::float_});
    check_meta(result.meta<1>(), {column_type::double_});
    BOOST_TEST(result.meta<2>().empty());

    // OK packet data
    BOOST_TEST(result.affected_rows<0>() == 1u);
    BOOST_TEST(result.last_insert_id<0>() == 2u);
    BOOST_TEST(result.warning_count<0>() == 3u);
    BOOST_TEST(result.info<0>() == "1st");
    BOOST_TEST(result.affected_rows<1>() == 4u);
    BOOST_TEST(result.last_insert_id<1>() == 5u);
    BOOST_TEST(result.warning_count<1>() == 6u);
    BOOST_TEST(result.info<1>() == "2nd");
    BOOST_TEST(result.affected_rows<2>() == 0u);
    BOOST_TEST(result.last_insert_id<2>() == 0u);
    BOOST_TEST(result.warning_count<2>() == 0u);
    BOOST_TEST(result.info<2>() == "3rd");
}

// Verify view validity
BOOST_AUTO_TEST_CASE(move_constructor)
{
    // Having this in heap makes easier to spot lifetime problems
    auto result = create_heap_results();

    // Obtain references
    auto rws0 = result->rows<0>();
    auto rws1 = result->rows<1>();
    auto meta0 = result->meta<0>();
    auto meta1 = result->meta<1>();
    auto info0 = result->info<0>();
    auto info1 = result->info<1>();

    // Move construct
    results_t result2(std::move(*result));
    result.reset();

    // Make sure that views are still valid
    BOOST_TEST((rws0[0] == row1{"abc", 10.1f}));
    BOOST_TEST((rws0[1] == row1{"def", 20.0f}));
    BOOST_TEST((rws1[0] == row2{42.0}));
    check_meta(meta0, {column_type::varchar, column_type::float_});
    check_meta(meta1, {column_type::double_});
    BOOST_TEST(info0 == "1st");
    BOOST_TEST(info1 == "2nd");

    // The new object holds the same data
    BOOST_TEST_REQUIRE(result2.has_value());
    BOOST_TEST((result2.rows<0>()[0] == row1{"abc", 10.1f}));
    check_meta(result2.meta<0>(), {column_type::varchar, column_type::float_});
    BOOST_TEST(result2.info<0>() == "1st");
}

BOOST_AUTO_TEST_CASE(move_assignment)
{
    // Having this in heap makes easier to spot lifetime problems
    auto result = create_heap_results();

    // Obtain references
    auto rws0 = result->rows<0>();
    auto rws1 = result->rows<1>();
    auto meta0 = result->meta<0>();
    auto meta1 = result->meta<1>();
    auto info0 = result->info<0>();
    auto info1 = result->info<1>();

    // Move construct
    results_t result2;
    result2 = std::move(*result);
    result.reset();

    // Make sure that views are still valid
    BOOST_TEST((rws0[0] == row1{"abc", 10.1f}));
    BOOST_TEST((rws0[1] == row1{"def", 20.0f}));
    BOOST_TEST((rws1[0] == row2{42.0}));
    check_meta(meta0, {column_type::varchar, column_type::float_});
    check_meta(meta1, {column_type::double_});
    BOOST_TEST(info0 == "1st");
    BOOST_TEST(info1 == "2nd");

    // The new object holds the same data
    BOOST_TEST_REQUIRE(result2.has_value());
    BOOST_TEST((result2.rows<0>()[0] == row1{"abc", 10.1f}));
    check_meta(result2.meta<0>(), {column_type::varchar, column_type::float_});
    BOOST_TEST(result2.info<0>() == "1st");
}

BOOST_AUTO_TEST_SUITE_END()

#endif
