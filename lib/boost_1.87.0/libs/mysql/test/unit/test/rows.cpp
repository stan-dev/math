//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/blob.hpp>
#include <boost/mysql/blob_view.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/rows.hpp>
#include <boost/mysql/rows_view.hpp>

#include <boost/test/unit_test.hpp>

#include <stdexcept>

#include "test_common/create_basic.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

BOOST_TEST_DONT_PRINT_LOG_VALUE(boost::mysql::detail::rows_iterator)

BOOST_AUTO_TEST_SUITE(test_rows)

// Check that references, pointers and iterators to certain rows' contents survive an operation
struct reference_checker
{
    rows::const_iterator it;
    rows_view rv;

    reference_checker(const rows& r) : it(r.begin()), rv(r) {}

    void check(const rows& new_row)
    {
        BOOST_TEST(it == new_row.begin());
        BOOST_TEST(rv == new_row);
    }
};

struct reference_checker_strs : reference_checker
{
    // indices in the matrix where a string/blob field resides
    std::size_t string_row_index, string_col_index, blob_row_index, blob_col_index;
    const char* string_ptr;
    const unsigned char* blob_ptr;

    static void assert_ptrs_equal(const char* ptr1, const char* ptr2)
    {
        // Otherwise, UTF thinks it's a C string, tries to print it and causes Valgrind errors
        BOOST_TEST(static_cast<const void*>(ptr1) == static_cast<const void*>(ptr2));
    }

    reference_checker_strs(
        const rows& r,
        std::size_t string_row_index,
        std::size_t string_col_index,
        std::size_t blob_row_index,
        std::size_t blob_col_index
    )
        : reference_checker(r),
          string_row_index(string_row_index),
          string_col_index(string_col_index),
          blob_row_index(blob_row_index),
          blob_col_index(blob_col_index),
          string_ptr(r.at(string_row_index).at(string_col_index).as_string().data()),
          blob_ptr(r.at(blob_row_index).at(blob_col_index).as_blob().data())
    {
    }

    void check(const rows& new_row)
    {
        reference_checker::check(new_row);
        assert_ptrs_equal(string_ptr, new_row.at(string_row_index).at(string_col_index).as_string().data());
        BOOST_TEST(blob_ptr == new_row.at(blob_row_index).at(blob_col_index).as_blob().data());
    }
};

BOOST_AUTO_TEST_CASE(default_ctor)
{
    rows r;
    BOOST_TEST(r.empty());
}

BOOST_AUTO_TEST_SUITE(ctor_from_rows_view)
BOOST_AUTO_TEST_CASE(empty)
{
    rows_view v;
    rows r(v);
    BOOST_TEST(r.empty());
}

BOOST_AUTO_TEST_CASE(non_empty)
{
    std::string s1("abc"), s2("");
    blob b{0x64, 0x10, 0x01};
    auto fields = make_fv_arr(s1, 1.0f, nullptr, s2, -1, b);
    auto v = makerowsv(fields.data(), fields.size(), 3);
    rows r(v);

    // r should be independent of the original fields/strings
    fields = make_fv_arr(0, 0, 0, 0, 0, 0);
    s1 = "other";
    s2 = "yet_another";
    b = {0x11, 0x12, 0x13};

    BOOST_TEST(r.size() == 2u);
    BOOST_TEST(r[0] == makerow("abc", 1.0f, nullptr));
    BOOST_TEST(r[1] == makerow("", -1, makebv("\x64\x10\1")));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(copy_ctor)
BOOST_AUTO_TEST_CASE(empty)
{
    rows r1;
    rows r2(r1);
    BOOST_TEST(r2.empty());
}

BOOST_AUTO_TEST_CASE(non_empty)
{
    rows r1 = makerows(3, "abc", 21.0f, makebv(""), "cdefg", 22.0f, makebv("\1\3\5"));
    rows r2(r1);
    r1 = makerows(2, 0, 0, 0, 0);  // r2 should be independent of r1

    BOOST_TEST(r2.size() == 2u);
    BOOST_TEST(r2[0] == makerow("abc", 21.0f, makebv("")));
    BOOST_TEST(r2[1] == makerow("cdefg", 22.0f, makebv("\1\3\5")));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(move_ctor)
BOOST_AUTO_TEST_CASE(empty)
{
    rows r1;

    // References, pointers, etc. should remain valid
    reference_checker refcheck(r1);

    rows r2(std::move(r1));
    r1 = makerows(1, "abc", 22);  // r2 should be independent of r1

    BOOST_TEST(r2.empty());
    refcheck.check(r2);
}

BOOST_AUTO_TEST_CASE(non_empty)
{
    rows r1 = makerows(3, makebv("abc"), 21.0f, "", makebv("\0\1\0"), 22.0f, "aaa");

    // References, pointers, etc. should remain valid
    reference_checker_strs refcheck(r1, 1, 2, 1, 0);

    rows r2(std::move(r1));
    r1 = makerows(2, 0, 0, 0, 0);  // r2 should be independent of r1

    BOOST_TEST(r2.size() == 2u);
    BOOST_TEST(r2[0] == makerow(makebv("abc"), 21.0f, ""));
    BOOST_TEST(r2[1] == makerow(makebv("\0\1\0"), 22.0f, "aaa"));
    refcheck.check(r2);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(copy_assignment)
BOOST_AUTO_TEST_CASE(empty)
{
    rows r1 = makerows(2, 42, "abcdef");
    rows r2;
    r1 = r2;
    r2 = makerows(2, 90, nullptr);  // r1 is independent of r2
    BOOST_TEST(r1.empty());
}

BOOST_AUTO_TEST_CASE(non_empty)
{
    rows r1 = makerows(1, 42, "abcdef", makebv("\0\2\1"));
    rows r2 = makerows(2, "a_very_long_string", nullptr, "", makebv("\7\0"), "cde", blob_view());
    r1 = r2;
    r2 = makerows(1, "another_string", 90, "yet_another");  // r1 is independent of r2

    BOOST_TEST(r1.size() == 3u);
    BOOST_TEST(r1[0] == makerow("a_very_long_string", nullptr));
    BOOST_TEST(r1[1] == makerow("", makebv("\7\0")));
    BOOST_TEST(r1[2] == makerow("cde", blob_view()));
}

BOOST_AUTO_TEST_CASE(self_assignment)
{
    rows r = makerows(2, "abc", 50u, "", makebv("abc"));
    const rows& ref = r;
    r = ref;

    BOOST_TEST(r.size() == 2u);
    BOOST_TEST(r[0] == makerow("abc", 50u));
    BOOST_TEST(r[1] == makerow("", makebv("abc")));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(move_assignment)
BOOST_AUTO_TEST_CASE(empty)
{
    rows r1 = makerows(1, 42, "abcdef");
    rows r2;

    // References, pointers, etc. should remain valid
    reference_checker refcheck(r2);
    r1 = std::move(r2);
    r2 = makerows(2, 90, nullptr);  // r1 is independent of r2

    BOOST_TEST(r1.empty());
    refcheck.check(r1);
}

BOOST_AUTO_TEST_CASE(non_empty)
{
    rows r1 = makerows(1, 42, "abcdef");
    rows r2 = makerows(3, "a_very_long_string", blob_view(), 50, "", makebv("\2\5\1"), 42);

    // References, pointers, etc. should remain valid
    reference_checker_strs refcheck(r2, 0, 0, 1, 1);

    r1 = std::move(r2);
    r2 = makerows(1, "another_string", 90, "yet_another");  // r1 is independent of r2

    BOOST_TEST(r1.size() == 2u);
    BOOST_TEST(r1[0] == makerow("a_very_long_string", blob_view(), 50));
    BOOST_TEST(r1[1] == makerow("", makebv("\2\5\1"), 42));
    refcheck.check(r1);
}

BOOST_AUTO_TEST_CASE(self_assignment)
{
    rows r = makerows(3, "abc", 50u, "fgh");
    rows&& ref = std::move(r);
    r = std::move(ref);  // this should leave r in a valid but unspecified state

    // r is in a valid but unspecified state; can be assigned to
    r = makerows(1, "abcdef");
    BOOST_TEST(r.size() == 1u);
    BOOST_TEST(r[0] == makerow("abcdef"));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(assignment_from_view)
BOOST_AUTO_TEST_CASE(empty)
{
    rows r = makerows(1, 42, "abcdef");
    r = rows_view();
    BOOST_TEST(r.empty());
    BOOST_TEST(r.num_columns() == 0u);
}

BOOST_AUTO_TEST_CASE(empty_different_num_columns)
{
    rows r;
    r = makerowsv(nullptr, 0, 2);

    BOOST_TEST(r.empty());
    BOOST_TEST(r.size() == 0u);
    BOOST_TEST(r.num_columns() == 2u);

    r = makerowsv(nullptr, 0, 3);

    BOOST_TEST(r.empty());
    BOOST_TEST(r.size() == 0u);
    BOOST_TEST(r.num_columns() == 3u);
}

BOOST_AUTO_TEST_CASE(non_empty)
{
    std::string s1("a_very_long_string"), s2("");
    blob b{0x78, 0x01, 0xff};
    rows r = makerows(1, 42, "abcdef", 90, "hij");
    auto fields = make_fv_arr(s1, nullptr, blob_view(), s2, "bec", b);
    r = makerowsv(fields.data(), fields.size(), 3);
    fields = make_fv_arr("abc", 42u, 9, 0, 0, 0);  // r should be independent of the original fields
    s1 = "another_string";                         // r should be independent of the original strings
    s2 = "yet_another";
    b = {0x00, 0xab, 0xcc, 0x78};

    BOOST_TEST_REQUIRE(r.size() == 2u);
    BOOST_TEST(r[0] == makerow("a_very_long_string", nullptr, blob_view()));
    BOOST_TEST(r[1] == makerow("", "bec", makebv("\x78\1\xff")));
    BOOST_TEST(r.num_columns() == 3u);
}

BOOST_AUTO_TEST_CASE(self_assignment)
{
    rows r = makerows(2, "abcdef", 42, "plk", "uv");
    r = rows_view(r);

    BOOST_TEST_REQUIRE(r.size() == 2u);
    BOOST_TEST(r[0] == makerow("abcdef", 42));
    BOOST_TEST(r[1] == makerow("plk", "uv"));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(at)
BOOST_AUTO_TEST_CASE(empty)
{
    rows v;
    BOOST_CHECK_THROW(v.at(0), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(one_column_one_row)
{
    rows r = makerows(1, 42u);
    BOOST_TEST(r.at(0) == makerow(42u));
    BOOST_CHECK_THROW(r.at(1), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(one_column_several_rows)
{
    rows r = makerows(1, 42u, "abc");
    BOOST_TEST(r.at(0) == makerow(42u));
    BOOST_TEST(r.at(1) == makerow("abc"));
    BOOST_CHECK_THROW(r.at(2), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(several_columns_one_row)
{
    rows r = makerows(2, 42u, "abc");
    BOOST_TEST(r.at(0) == makerow(42u, "abc"));
    BOOST_CHECK_THROW(r.at(1), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(several_columns_several_rows)
{
    rows r = makerows(2, 42u, "abc", nullptr, "bcd", 90u, nullptr);
    BOOST_TEST(r.at(0) == makerow(42u, "abc"));
    BOOST_TEST(r.at(1) == makerow(nullptr, "bcd"));
    BOOST_TEST(r.at(2) == makerow(90u, nullptr));
    BOOST_CHECK_THROW(r.at(3), std::out_of_range);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(operator_square_brackets)
BOOST_AUTO_TEST_CASE(one_column_one_row)
{
    rows r = makerows(1, 42u);
    BOOST_TEST(r[0] == makerow(42u));
}

BOOST_AUTO_TEST_CASE(one_column_several_rows)
{
    rows r = makerows(1, 42u, "abc");
    BOOST_TEST(r[0] == makerow(42u));
    BOOST_TEST(r[1] == makerow("abc"));
}

BOOST_AUTO_TEST_CASE(several_columns_one_row)
{
    rows r = makerows(2, 42u, "abc");
    BOOST_TEST(r[0] == makerow(42u, "abc"));
}

BOOST_AUTO_TEST_CASE(several_columns_several_rows)
{
    rows r = makerows(2, 42u, "abc", nullptr, "bcd", 90u, nullptr);
    BOOST_TEST(r[0] == makerow(42u, "abc"));
    BOOST_TEST(r[1] == makerow(nullptr, "bcd"));
    BOOST_TEST(r[2] == makerow(90u, nullptr));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_CASE(front)
{
    rows r = makerows(2, 42u, "abc", nullptr, "bcde");
    BOOST_TEST(r.front() == makerow(42u, "abc"));
}

BOOST_AUTO_TEST_CASE(back)
{
    rows r = makerows(2, 70.0f, "abc", nullptr, "bcde");
    BOOST_TEST(r.back() == makerow(nullptr, "bcde"));
}

BOOST_AUTO_TEST_CASE(empty)
{
    BOOST_TEST(rows().empty());
    BOOST_TEST(!makerows(1, 42u).empty());
}

BOOST_AUTO_TEST_SUITE(size)
BOOST_AUTO_TEST_CASE(empty)
{
    rows r;
    BOOST_TEST(r.size() == 0u);
}

BOOST_AUTO_TEST_CASE(one_column_one_row)
{
    rows r = makerows(1, 42u);
    BOOST_TEST(r.size() == 1u);
}

BOOST_AUTO_TEST_CASE(one_column_several_rows)
{
    rows r = makerows(1, 42u, "abc");
    BOOST_TEST(r.size() == 2u);
}

BOOST_AUTO_TEST_CASE(several_columns_one_row)
{
    rows r = makerows(2, 42u, "abc");
    BOOST_TEST(r.size() == 1u);
}

BOOST_AUTO_TEST_CASE(several_columns_several_rows)
{
    rows r = makerows(3, 42u, "abc", nullptr, "bcd", 90u, nullptr);
    BOOST_TEST(r.size() == 2u);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(operator_rows_view)
BOOST_AUTO_TEST_CASE(empty)
{
    rows r;
    auto rv = static_cast<rows_view>(r);
    BOOST_TEST(rv.size() == 0u);
}

BOOST_AUTO_TEST_CASE(non_empty)
{
    rows r = makerows(3, 42u, 4.2f, "abcde", 90u, nullptr, "def");
    auto rv = static_cast<rows_view>(r);
    BOOST_TEST(rv.size() == 2u);
    BOOST_TEST(rv[0] == makerow(42u, 4.2f, "abcde"));
    BOOST_TEST(rv[1] == makerow(90u, nullptr, "def"));
}

BOOST_AUTO_TEST_CASE(cleared)
{
    rows r = makerows(3, 42u, 4.2f, "abcde", 90u, nullptr, "def");
    r = rows();
    auto rv = static_cast<rows_view>(r);
    BOOST_TEST(rv.empty());
    BOOST_TEST(rv.size() == 0u);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
