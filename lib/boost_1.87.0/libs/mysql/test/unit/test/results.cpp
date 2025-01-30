//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/column_type.hpp>
#include <boost/mysql/results.hpp>

#include <boost/test/unit_test.hpp>

#include <stdexcept>

#include "test_common/check_meta.hpp"
#include "test_common/create_basic.hpp"
#include "test_unit/create_execution_processor.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

BOOST_TEST_DONT_PRINT_LOG_VALUE(results::iterator)

BOOST_AUTO_TEST_SUITE(test_results)

results create_initial_results()
{
    results res;
    exec_access(get_iface(res))
        .meta({meta_builder().type(column_type::varchar).build_coldef()})
        .row("abc")
        .row(nullptr)
        .ok(ok_builder().affected_rows(1).last_insert_id(2).warnings(3).info("1st").more_results(true).build()
        )
        .meta({meta_builder().type(column_type::tinyint).build_coldef()})
        .row(42)
        .ok(ok_builder()
                .affected_rows(4)
                .last_insert_id(5)
                .warnings(6)
                .info("2nd")
                .out_params(true)
                .more_results(true)
                .build())
        .ok(ok_builder().info("3rd").build());
    return res;
}

struct fixture
{
    results result{create_initial_results()};
};

BOOST_AUTO_TEST_CASE(has_value)
{
    // Default construction
    results result;
    BOOST_TEST_REQUIRE(!result.has_value());

    // With value
    result = create_initial_results();
    BOOST_TEST_REQUIRE(result.has_value());
}

BOOST_AUTO_TEST_SUITE(iterators)

BOOST_FIXTURE_TEST_CASE(basic, fixture)
{
    // Obtain iterators
    auto it = result.begin();   // should point to resultset 0
    auto itend = result.end();  // should point to resultset 3 (1 past end)

    // Check dereference
    BOOST_TEST((*it).info() == "1st");
    BOOST_TEST(it->info() == "1st");

    // Check ==
    BOOST_TEST(!(it == itend));
    BOOST_TEST(!(itend == it));
    BOOST_TEST(it == result.begin());
    BOOST_TEST(it == it);
    BOOST_TEST(itend == result.end());
    BOOST_TEST(itend == itend);

    // Check !=
    BOOST_TEST(it != itend);
    BOOST_TEST(itend != it);
    BOOST_TEST(!(it != result.begin()));
    BOOST_TEST(!(it != it));
    BOOST_TEST(!(itend != result.end()));
    BOOST_TEST(!(itend != itend));
}

BOOST_FIXTURE_TEST_CASE(prefix_increment, fixture)
{
    auto it = result.begin();
    auto& ref = (++it);
    BOOST_TEST(&ref == &it);
    BOOST_TEST(it->info() == "2nd");
    BOOST_TEST(it == result.begin() + 1);
}

BOOST_FIXTURE_TEST_CASE(postfix_increment, fixture)
{
    auto it = result.begin();
    auto it2 = it++;
    BOOST_TEST(it2 == result.begin());
    BOOST_TEST(it == result.begin() + 1);
    BOOST_TEST(it->info() == "2nd");
}

BOOST_FIXTURE_TEST_CASE(prefix_decrement, fixture)
{
    auto it = result.end();
    auto& ref = (--it);
    BOOST_TEST(&ref == &it);
    BOOST_TEST(it->info() == "3rd");
    BOOST_TEST(it == result.begin() + 2);
}

BOOST_FIXTURE_TEST_CASE(postfix_decrement, fixture)
{
    auto it = result.end();
    auto it2 = it--;
    BOOST_TEST(it2 == result.end());
    BOOST_TEST(it == result.begin() + 2);
    BOOST_TEST(it->info() == "3rd");
}

BOOST_FIXTURE_TEST_CASE(operator_square_brackets, fixture)
{
    auto it = result.begin();
    BOOST_TEST(it[0].info() == "1st");
    BOOST_TEST(it[1].info() == "2nd");
    BOOST_TEST(it[2].info() == "3rd");
}

BOOST_FIXTURE_TEST_CASE(operator_plus, fixture)
{
    auto it = result.begin();

    // Increment by 1
    auto it2 = it + 1;
    BOOST_TEST(it2->info() == "2nd");

    // Reversed operands
    it2 = 1 + it2;
    BOOST_TEST(it2->info() == "3rd");

    // Increment by more than 1
    BOOST_TEST(result.begin() + 3 == result.end());

    // Increment by 0
    BOOST_TEST(result.begin() + 0 == result.begin());

    // Negative increment
    BOOST_TEST(result.end() + (-2) == result.begin() + 1);
}

BOOST_FIXTURE_TEST_CASE(operator_plus_equals, fixture)
{
    auto it = result.begin();

    // Increment by 1
    it += 1;
    BOOST_TEST(it->info() == "2nd");

    // Increment by more than
    it += 2;
    BOOST_TEST(it == result.end());

    // Increment by 0
    it += 0;
    BOOST_TEST(it == result.end());

    // Negative increment
    it += (-2);
    BOOST_TEST(it == result.begin() + 1);
}

BOOST_FIXTURE_TEST_CASE(operator_minus, fixture)
{
    auto it = result.end();

    // Decrement by 1
    auto it2 = it - 1;
    BOOST_TEST(it2->info() == "3rd");

    // Decrement by more than 1
    BOOST_TEST(result.end() - 3 == result.begin());

    // Decrement by 0
    BOOST_TEST(result.end() - 0 == result.end());

    // Negative decrement
    BOOST_TEST(result.begin() - (-2) == result.begin() + 2);
}

BOOST_FIXTURE_TEST_CASE(operator_minus_equals, fixture)
{
    auto it = result.end();

    // Decrement by 1
    it -= 1;
    BOOST_TEST(it->info() == "3rd");

    // Decrement by more than 1
    it -= 2;
    BOOST_TEST(it == result.begin());

    // Decrement by 0
    it -= 0;
    BOOST_TEST(it == result.begin());

    // Negative decrement
    it -= (-2);
    BOOST_TEST(it == result.begin() + 2);
}

BOOST_FIXTURE_TEST_CASE(difference, fixture)
{
    auto first = result.begin();
    auto second = result.begin() + 1;
    auto last = result.end();

    BOOST_TEST(last - first == 3);
    BOOST_TEST(last - second == 2);
    BOOST_TEST(last - last == 0);
    BOOST_TEST(first - last == -3);
    BOOST_TEST(second - last == -2);
    BOOST_TEST(last - last == 0);
    BOOST_TEST(first - first == 0);
}

BOOST_FIXTURE_TEST_CASE(relational, fixture)
{
    auto first = result.begin();
    auto second = result.begin() + 1;
    auto third = result.begin() + 2;

    // Less than
    BOOST_TEST(first < second);
    BOOST_TEST(first <= second);
    BOOST_TEST(!(first > second));
    BOOST_TEST(!(first >= second));

    // Equal
    BOOST_TEST(!(second < second));
    BOOST_TEST(second <= second);
    BOOST_TEST(!(second > second));
    BOOST_TEST(second >= second);

    // Greater than
    BOOST_TEST(!(third < second));
    BOOST_TEST(!(third <= second));
    BOOST_TEST(third > second);
    BOOST_TEST(third >= second);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_CASE(collection_fns, fixture)
{
    // at
    BOOST_TEST(result.at(0).info() == "1st");
    BOOST_TEST(result.at(1).info() == "2nd");
    BOOST_TEST(result.at(2).info() == "3rd");
    BOOST_CHECK_THROW(result.at(3), std::out_of_range);

    // operator[]
    BOOST_TEST(result[0].info() == "1st");
    BOOST_TEST(result[1].info() == "2nd");
    BOOST_TEST(result[2].info() == "3rd");

    // front & back
    BOOST_TEST(result.front().info() == "1st");
    BOOST_TEST(result.back().info() == "3rd");

    // size & empty
    BOOST_TEST(result.size() == 3u);
    BOOST_TEST(!result.empty());
}

// Verify view validity
BOOST_AUTO_TEST_CASE(move_constructor)
{
    // Having this in heap helps spot lifetime issues
    std::unique_ptr<results> result{new results(create_initial_results())};

    // Obtain references. Note that iterators and resultset_view's don't remain valid.
    auto rws = result->rows();
    auto meta = result->meta();
    auto info = result->info();

    // Move construct
    results result2(std::move(*result));
    result.reset();

    // Make sure that views are still valid
    BOOST_TEST(rws == makerows(1, "abc", nullptr));
    check_meta(meta, {column_type::varchar});
    BOOST_TEST(info == "1st");

    // The new object holds the same data
    BOOST_TEST_REQUIRE(result2.has_value());
    BOOST_TEST(result2.rows() == makerows(1, "abc", nullptr));
    check_meta(result2.meta(), {column_type::varchar});
    BOOST_TEST(result2.info() == "1st");
}

BOOST_AUTO_TEST_CASE(move_assignment)
{
    // Having this in heap helps spot lifetime issues
    std::unique_ptr<results> result{new results(create_initial_results())};

    // Obtain references
    auto rws = result->rows();
    auto meta = result->meta();
    auto info = result->info();

    // Move construct
    results result2;
    result2 = std::move(*result);
    result.reset();

    // Make sure that views are still valid
    BOOST_TEST(rws == makerows(1, "abc", nullptr));
    check_meta(meta, {column_type::varchar});
    BOOST_TEST(info == "1st");

    // The new object holds the same data
    BOOST_TEST_REQUIRE(result2.has_value());
    BOOST_TEST(result2.rows() == makerows(1, "abc", nullptr));
    check_meta(result2.meta(), {column_type::varchar});
    BOOST_TEST(result2.info() == "1st");
}

BOOST_AUTO_TEST_SUITE_END()
