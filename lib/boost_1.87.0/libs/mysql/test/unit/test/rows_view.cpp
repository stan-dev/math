//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/field_view.hpp>
#include <boost/mysql/rows_view.hpp>

#include <boost/test/unit_test.hpp>

#include <stdexcept>

#include "test_common/create_basic.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

BOOST_AUTO_TEST_SUITE(test_rows_view)

BOOST_AUTO_TEST_CASE(default_ctor)
{
    rows_view v;
    BOOST_TEST(v.empty());
}

// Regression checks. The cases below caused segfaults in debug builds
BOOST_AUTO_TEST_SUITE(init_ctor)
BOOST_AUTO_TEST_CASE(fieldsnonnull_nfieldszero_ncolsnonzero)
{
    field_view fv;
    auto v = makerowsv(&fv, 0, 3);
    BOOST_TEST(v.empty());
    BOOST_TEST(v.size() == 0u);
    BOOST_TEST(v.num_columns() == 3u);
    BOOST_CHECK_THROW(v.at(0), std::out_of_range);
    std::vector<row_view> vec(v.begin(), v.end());
    BOOST_TEST(vec.empty());
}

BOOST_AUTO_TEST_CASE(fieldsnonnull_nfieldszero_ncolszero)
{
    field_view fv;
    auto v = makerowsv(&fv, 0, 0);
    BOOST_TEST(v.empty());
    BOOST_TEST(v.size() == 0u);
    BOOST_TEST(v.num_columns() == 0u);
    BOOST_CHECK_THROW(v.at(0), std::out_of_range);
    std::vector<row_view> vec(v.begin(), v.end());
    BOOST_TEST(vec.empty());
}

BOOST_AUTO_TEST_CASE(fieldsnull_nfieldszero_ncolsnonzero)
{
    auto v = makerowsv(nullptr, 0, 2);
    BOOST_TEST(v.empty());
    BOOST_TEST(v.size() == 0u);
    BOOST_TEST(v.num_columns() == 2u);
    BOOST_CHECK_THROW(v.at(0), std::out_of_range);
    std::vector<row_view> vec(v.begin(), v.end());
    BOOST_TEST(vec.empty());
}

BOOST_AUTO_TEST_CASE(fieldsnull_nfieldszero_ncolszero)
{
    auto v = makerowsv(nullptr, 0, 0);
    BOOST_TEST(v.empty());
    BOOST_TEST(v.size() == 0u);
    BOOST_TEST(v.num_columns() == 0u);
    BOOST_CHECK_THROW(v.at(0), std::out_of_range);
    std::vector<row_view> vec(v.begin(), v.end());
    BOOST_TEST(vec.empty());
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(at)
BOOST_AUTO_TEST_CASE(empty)
{
    rows_view v;
    BOOST_CHECK_THROW(v.at(0), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(one_column_one_row)
{
    auto fields = make_fv_arr(42u);
    auto v = makerowsv(fields.data(), 1, 1);
    BOOST_TEST(v.at(0) == makerow(42u));
    BOOST_CHECK_THROW(v.at(1), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(one_column_several_rows)
{
    auto fields = make_fv_arr(42u, "abc");
    auto v = makerowsv(fields.data(), 2, 1);
    BOOST_TEST(v.at(0) == makerow(42u));
    BOOST_TEST(v.at(1) == makerow("abc"));
    BOOST_CHECK_THROW(v.at(2), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(several_columns_one_row)
{
    auto fields = make_fv_arr(42u, "abc");
    auto v = makerowsv(fields.data(), 2, 2);
    BOOST_TEST(v.at(0) == makerow(42u, "abc"));
    BOOST_CHECK_THROW(v.at(1), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(several_columns_several_rows)
{
    auto fields = make_fv_arr(42u, "abc", nullptr, "bcd", 90u, nullptr);
    auto v = makerowsv(fields.data(), 6, 2);
    BOOST_TEST(v.at(0) == makerow(42u, "abc"));
    BOOST_TEST(v.at(1) == makerow(nullptr, "bcd"));
    BOOST_TEST(v.at(2) == makerow(90u, nullptr));
    BOOST_CHECK_THROW(v.at(3), std::out_of_range);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(operator_square_brackets)
BOOST_AUTO_TEST_CASE(one_column_one_row)
{
    auto fields = make_fv_arr(42u);
    auto v = makerowsv(fields.data(), 1, 1);
    BOOST_TEST(v[0] == makerow(42u));
}

BOOST_AUTO_TEST_CASE(one_column_several_rows)
{
    auto fields = make_fv_arr(42u, "abc");
    auto v = makerowsv(fields.data(), 2, 1);
    BOOST_TEST(v[0] == makerow(42u));
    BOOST_TEST(v[1] == makerow("abc"));
}

BOOST_AUTO_TEST_CASE(several_columns_one_row)
{
    auto fields = make_fv_arr(42u, "abc");
    auto v = makerowsv(fields.data(), 2, 2);
    BOOST_TEST(v[0] == makerow(42u, "abc"));
}

BOOST_AUTO_TEST_CASE(several_columns_several_rows)
{
    auto fields = make_fv_arr(42u, "abc", nullptr, "bcd", 90u, nullptr);
    auto v = makerowsv(fields.data(), 6, 2);
    BOOST_TEST(v[0] == makerow(42u, "abc"));
    BOOST_TEST(v[1] == makerow(nullptr, "bcd"));
    BOOST_TEST(v[2] == makerow(90u, nullptr));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_CASE(front)
{
    auto fields = make_fv_arr(42u, "abc", nullptr, "bcde");
    auto v = makerowsv(fields.data(), 4, 2);
    BOOST_TEST(v.front() == makerow(42u, "abc"));
}

BOOST_AUTO_TEST_CASE(back)
{
    auto fields = make_fv_arr(42u, "abc", nullptr, "bcde");
    auto v = makerowsv(fields.data(), 4, 2);
    BOOST_TEST(v.back() == makerow(nullptr, "bcde"));
}

BOOST_AUTO_TEST_CASE(empty)
{
    BOOST_TEST(rows_view().empty());

    auto fields = make_fv_arr(42u);
    BOOST_TEST(!makerowsv(fields.data(), 1, 1).empty());
}

BOOST_AUTO_TEST_SUITE(size)
BOOST_AUTO_TEST_CASE(empty)
{
    rows_view v;
    BOOST_TEST(v.size() == 0u);
}

BOOST_AUTO_TEST_CASE(one_column_one_row)
{
    auto fields = make_fv_arr(42u);
    auto v = makerowsv(fields.data(), 1, 1);
    BOOST_TEST(v.size() == 1u);
}

BOOST_AUTO_TEST_CASE(one_column_several_rows)
{
    auto fields = make_fv_arr(42u, "abc");
    auto v = makerowsv(fields.data(), 2, 1);
    BOOST_TEST(v.size() == 2u);
}

BOOST_AUTO_TEST_CASE(several_columns_one_row)
{
    auto fields = make_fv_arr(42u, "abc");
    auto v = makerowsv(fields.data(), 2, 2);
    BOOST_TEST(v.size() == 1u);
}

BOOST_AUTO_TEST_CASE(several_columns_several_rows)
{
    auto fields = make_fv_arr(42u, "abc", nullptr, "bcd", 90u, nullptr);
    auto v = makerowsv(fields.data(), 6, 2);
    BOOST_TEST(v.size() == 3u);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
