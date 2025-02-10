//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/field_view.hpp>
#include <boost/mysql/statement.hpp>

#include <boost/test/unit_test.hpp>

#include <tuple>

#include "test_common/create_basic.hpp"
#include "test_unit/create_statement.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

BOOST_AUTO_TEST_SUITE(test_statement)

BOOST_AUTO_TEST_CASE(default_ctor)
{
    statement stmt;
    BOOST_TEST(!stmt.valid());
}

BOOST_AUTO_TEST_CASE(member_fns)
{
    auto stmt = statement_builder().id(1).num_params(3).build();

    BOOST_TEST(stmt.valid());
    BOOST_TEST(stmt.num_params() == 3u);
    BOOST_TEST(stmt.id() == 1u);
}

statement create_valid_stmt() { return statement_builder().id(1).build(); }

BOOST_AUTO_TEST_SUITE(bind_field_like)
BOOST_AUTO_TEST_CASE(regular)
{
    std::string s("def");
    blob blb;
    auto b = create_valid_stmt().bind(42, std::string("abc"), std::ref(s), s, blb);
    using tup_type = std::tuple<int, std::string, std::string&, std::string, blob>;
    static_assert(std::is_same<decltype(b), bound_statement_tuple<tup_type>>::value, "");
}

BOOST_AUTO_TEST_CASE(empty)
{
    auto b = create_valid_stmt().bind();
    static_assert(std::is_same<decltype(b), bound_statement_tuple<std::tuple<>>>::value, "");
}

BOOST_AUTO_TEST_CASE(stmt_const)
{
    const statement stmt = create_valid_stmt();
    auto b = stmt.bind();  // compiles
    static_assert(std::is_same<decltype(b), bound_statement_tuple<std::tuple<>>>::value, "");
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(bind_field_like_tuple)
BOOST_AUTO_TEST_CASE(regular)
{
    auto b = create_valid_stmt().bind(std::make_tuple(42, 4.2f));
    using expected_type = bound_statement_tuple<std::tuple<int, float>>;
    static_assert(std::is_same<decltype(b), expected_type>::value, "");
}

BOOST_AUTO_TEST_CASE(tuple_const_reference)
{
    const auto params = std::make_tuple(42, 4.2f);
    auto b = create_valid_stmt().bind(params);
    using expected_type = bound_statement_tuple<std::tuple<int, float>>;
    static_assert(std::is_same<decltype(b), expected_type>::value, "");
}

BOOST_AUTO_TEST_CASE(tuple_reference)
{
    auto params = std::make_tuple(42, 4.2f);
    auto b = create_valid_stmt().bind(params);
    using expected_type = bound_statement_tuple<std::tuple<int, float>>;
    static_assert(std::is_same<decltype(b), expected_type>::value, "");
}

BOOST_AUTO_TEST_CASE(stmt_const)
{
    const statement stmt = create_valid_stmt();
    auto b = stmt.bind(std::make_tuple(42));  // compiles
    using expected_type = bound_statement_tuple<std::tuple<int>>;
    static_assert(std::is_same<decltype(b), expected_type>::value, "");
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(bind_iterator_range)
BOOST_AUTO_TEST_CASE(regular)
{
    auto fvs = make_fv_arr(42, "abc");
    auto b = create_valid_stmt().bind(fvs.begin(), fvs.end());
    using expected_type = bound_statement_iterator_range<std::array<field_view, 2>::iterator>;
    static_assert(std::is_same<decltype(b), expected_type>::value, "");
}

BOOST_AUTO_TEST_CASE(stmt_const)
{
    const statement stmt = create_valid_stmt();
    auto fvs = make_fv_arr(42, "abc");
    auto b = stmt.bind(fvs.begin(), fvs.end());  // compiles
    using expected_type = bound_statement_iterator_range<std::array<field_view, 2>::iterator>;
    static_assert(std::is_same<decltype(b), expected_type>::value, "");
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
