//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/blob.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/row.hpp>
#include <boost/mysql/row_view.hpp>

#include <boost/test/unit_test.hpp>

#include <string>
#include <vector>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/create_basic.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

BOOST_AUTO_TEST_SUITE(test_row)

// Check that references, pointers and iterators to certain row's contents survive an operation
struct reference_checker
{
    row::const_iterator it;
    row_view rv;

    reference_checker(const row& r) : it(r.begin()), rv(r) {}

    void check(const row& new_row)
    {
        BOOST_TEST(it == new_row.begin());
        BOOST_TEST(rv == new_row);
    }
};

struct reference_checker_strs : reference_checker
{
    std::size_t string_index, blob_index;  // indices in the row where a string/blob field resides
    const char* string_ptr;
    const unsigned char* blob_ptr;

    static void assert_ptrs_equal(const char* ptr1, const char* ptr2)
    {
        // Otherwise, UTF thinks it's a C string, tries to print it and causes Valgrind errors
        BOOST_TEST(static_cast<const void*>(ptr1) == static_cast<const void*>(ptr2));
    }

    reference_checker_strs(const row& r, std::size_t string_index, std::size_t blob_index)
        : reference_checker(r),
          string_index(string_index),
          blob_index(blob_index),
          string_ptr(r.at(string_index).as_string().data()),
          blob_ptr(r.at(blob_index).as_blob().data())
    {
    }

    void check(const row& new_row)
    {
        reference_checker::check(new_row);
        assert_ptrs_equal(string_ptr, new_row.at(string_index).as_string().data());
        BOOST_TEST(blob_ptr == new_row.at(blob_index).as_blob().data());
    }
};

BOOST_AUTO_TEST_CASE(default_ctor)
{
    row r;
    BOOST_TEST(r.empty());
}

BOOST_AUTO_TEST_SUITE(ctor_from_frow_view)
BOOST_AUTO_TEST_CASE(empty)
{
    row_view v;
    row r(v);
    BOOST_TEST(r.empty());
}

BOOST_AUTO_TEST_CASE(non_empty)
{
    std::string s1("test"), s2("");
    blob b{0x00, 0xab, 0xf5};
    auto fields = make_fv_arr(42, s1, 5.0f, b, s2);
    row r(makerowv(fields.data(), fields.size()));

    // Fields still valid even when the original source of the view changed
    fields = make_fv_arr(0, 0, 0, 0, 0);
    s1 = "other";
    s2 = "abcdef";
    b = {0xff, 0xa4, 0x02};

    BOOST_TEST(r.size() == 5u);
    BOOST_TEST(r[0] == field_view(42));
    BOOST_TEST(r[1] == field_view("test"));
    BOOST_TEST(r[2] == field_view(5.0f));
    BOOST_TEST(r[3] == field_view(makebv("\0\xab\xf5")));
    BOOST_TEST(r[4] == field_view(""));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(copy_ctor)
BOOST_AUTO_TEST_CASE(empty)
{
    row r1;
    row r2(r1);
    r1 = makerow(42, "test");  // r2 should be independent of r1

    BOOST_TEST(r2.empty());
}

BOOST_AUTO_TEST_CASE(non_empty)
{
    row r1 = makerow("", 42, "test", makebv("\0\3\2"));
    row r2(r1);
    r1 = makerow("another_string", 4.2f, "");  // r2 should be independent of r1

    BOOST_TEST(r2.size() == 4u);
    BOOST_TEST(r2[0] == field_view(""));
    BOOST_TEST(r2[1] == field_view(42));
    BOOST_TEST(r2[2] == field_view("test"));
    BOOST_TEST(r2[3] == field_view(makebv("\0\3\2")));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(move_ctor)
BOOST_AUTO_TEST_CASE(empty)
{
    row r1;

    // References, pointers, etc. should remain valid
    reference_checker refcheck(r1);

    row r2(std::move(r1));
    r1 = makerow(42, "test");  // r2 should be independent of r1

    BOOST_TEST(r2.empty());
    refcheck.check(r2);
}

BOOST_AUTO_TEST_CASE(non_empty)
{
    row r1 = makerow("", 42, "test", makebv("\0\5\xff"));

    // References, pointers, etc should remain valid
    reference_checker_strs refcheck(r1, 2, 3);

    // Move
    row r2(std::move(r1));
    r1 = makerow("another_string", 4.2f, "", makebv("\1\5\xab"));  // r2 should be independent of r1

    BOOST_TEST(r2.size() == 4u);
    BOOST_TEST(r2[0] == field_view(""));
    BOOST_TEST(r2[1] == field_view(42));
    BOOST_TEST(r2[2] == field_view("test"));
    BOOST_TEST(r2[3] == field_view(makebv("\0\5\xff")));

    // References, pointers, etc still valid
    refcheck.check(r2);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(copy_assignment)
BOOST_AUTO_TEST_CASE(empty)
{
    row r1 = makerow(42, "abcdef");
    row r2;
    r1 = r2;
    r2 = makerow(90, nullptr);  // r1 is independent of r2
    BOOST_TEST(r1.empty());
}

BOOST_AUTO_TEST_CASE(non_empty)
{
    row r1 = makerow(42, "abcdef", makebv("\0\1\2"));
    row r2 = makerow("a_very_long_string", nullptr, "", makebv("\3\4\5"));
    r1 = r2;
    r2 = makerow("another_string", 90, "yet_another");  // r1 is independent of r2

    BOOST_TEST(r1.size() == 4u);
    BOOST_TEST(r1[0] == field_view("a_very_long_string"));
    BOOST_TEST(r1[1] == field_view());
    BOOST_TEST(r1[2] == field_view(""));
    BOOST_TEST(r1[3] == field_view(makebv("\3\4\5")));
}

BOOST_AUTO_TEST_CASE(self_assignment)
{
    row r = makerow("abc", 50u, "fgh");
    const row& ref = r;
    r = ref;

    BOOST_TEST(r.size() == 3u);
    BOOST_TEST(r[0] == field_view("abc"));
    BOOST_TEST(r[1] == field_view(50u));
    BOOST_TEST(r[2] == field_view("fgh"));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(move_assignment)
BOOST_AUTO_TEST_CASE(empty)
{
    row r1 = makerow(42, "abcdef");
    row r2;
    row_view rv(r2);
    r1 = std::move(r2);
    r2 = makerow(90, nullptr);  // r1 is independent of r2
    BOOST_TEST(r1.empty());
    BOOST_TEST(rv == r1);
}

BOOST_AUTO_TEST_CASE(non_empty)
{
    row r1 = makerow(42, "abcdef", makebv("\0\4\1"));
    row r2 = makerow("a_very_long_string", nullptr, "", makebv("\7\1\2"));

    // References, pointers, etc should remain valid
    reference_checker_strs refcheck(r2, 0, 3);

    // Move
    r1 = std::move(r2);
    r2 = makerow("another_string", 90, "yet_another", makebv("\0\0"));  // r1 is independent of r2

    BOOST_TEST(r1.size() == 4u);
    BOOST_TEST(r1[0] == field_view("a_very_long_string"));
    BOOST_TEST(r1[1] == field_view());
    BOOST_TEST(r1[2] == field_view(""));
    BOOST_TEST(r1[3] == field_view(makebv("\7\1\2")));

    // References, pointers, etc still valid
    refcheck.check(r1);
}

BOOST_AUTO_TEST_CASE(self_assignment)
{
    row r = makerow("abc", 50u, "fgh", makebv("\0\4"));
    row&& ref = std::move(r);
    r = std::move(ref);  // this should leave r in a valid but unspecified state

    // r is in a valid but unspecified state; can be assigned to
    r = makerow("abcdef");
    BOOST_TEST(r.size() == 1u);
    BOOST_TEST(r[0] == field_view("abcdef"));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(assignment_from_view)
BOOST_AUTO_TEST_CASE(empty)
{
    row r = makerow(42, "abcdef", makebv("\0\xae"));
    r = row_view();
    BOOST_TEST(r.empty());
}

BOOST_AUTO_TEST_CASE(non_empty)
{
    std::string s1("a_very_long_string"), s2("");
    blob b{0x00, 0xfa};
    row r = makerow(42, "abcdef", makebv("\0\1"));
    auto fields = make_fv_arr(s1, nullptr, s2, b);
    r = makerowv(fields.data(), fields.size());
    fields = make_fv_arr("abc", 42u, 9, nullptr);  // r should be independent of the original fields
    s1 = "another_string";                         // r should be independent of the original strings
    s2 = "yet_another";
    b = {0xac, 0x32, 0x21, 0x50};

    BOOST_TEST(r.size() == 4u);
    BOOST_TEST(r[0] == field_view("a_very_long_string"));
    BOOST_TEST(r[1] == field_view());
    BOOST_TEST(r[2] == field_view(""));
    BOOST_TEST(r[3] == field_view(makebv("\0\xfa")));
}

BOOST_AUTO_TEST_CASE(self_assignment)
{
    row r = makerow("abcdef", 42, "plk");
    r = row_view(r);

    BOOST_TEST(r.size() == 3u);
    BOOST_TEST(r[0] == field_view("abcdef"));
    BOOST_TEST(r[1] == field_view(42));
    BOOST_TEST(r[2] == field_view("plk"));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(at)
BOOST_AUTO_TEST_CASE(empty)
{
    row r;
    BOOST_CHECK_THROW(r.at(0), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(in_range)
{
    row r = makerow(42, 50u, "test");
    BOOST_TEST(r.at(0) == field_view(42));
    BOOST_TEST(r.at(1) == field_view(50u));
    BOOST_TEST(r.at(2) == field_view("test"));
}

BOOST_AUTO_TEST_CASE(out_of_range)
{
    row r = makerow(42, 50u, "test");
    BOOST_CHECK_THROW(r.at(3), std::out_of_range);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_CASE(front)
{
    auto r = makerow(42, 50u, "test");
    BOOST_TEST(r.front() == field_view(42));
}

BOOST_AUTO_TEST_CASE(back)
{
    BOOST_TEST(makerow(42, 50u, "test").front() == field_view(42));
    BOOST_TEST(makerow(42).back() == field_view(42));
}

BOOST_AUTO_TEST_CASE(empty)
{
    BOOST_TEST(row().empty());
    BOOST_TEST(!makerow(42).empty());
    BOOST_TEST(!makerow(42, 50u).empty());
}

BOOST_AUTO_TEST_CASE(size)
{
    BOOST_TEST(row().size() == 0u);
    BOOST_TEST(makerow(42).size() == 1u);
    BOOST_TEST(makerow(50, nullptr).size() == 2u);
}

// As iterators are regular pointers, we don't perform
// exhaustive testing on iteration
BOOST_AUTO_TEST_SUITE(iterators)
BOOST_AUTO_TEST_CASE(empty)
{
    const row r{};  // can be called on const objects
    BOOST_TEST(r.begin() == nullptr);
    BOOST_TEST(r.end() == nullptr);
    std::vector<field_view> vec{r.begin(), r.end()};
    BOOST_TEST(vec.empty());
}

BOOST_AUTO_TEST_CASE(multiple_elms)
{
    const row r = makerow(42, 50u, "test");  // can be called on const objects
    BOOST_TEST(r.begin() != nullptr);
    BOOST_TEST(r.end() != nullptr);
    BOOST_TEST(std::distance(r.begin(), r.end()) == 3);

    std::vector<field_view> vec{r.begin(), r.end()};
    BOOST_TEST(vec.size() == 3u);
    BOOST_TEST(vec[0] == field_view(42));
    BOOST_TEST(vec[1] == field_view(50u));
    BOOST_TEST(vec[2] == field_view("test"));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(operator_row_view)
BOOST_AUTO_TEST_CASE(empty)
{
    row r;
    row_view rv(r);

    BOOST_TEST(rv.empty());
    BOOST_TEST(rv.size() == 0u);
    BOOST_TEST(rv.begin() == nullptr);
    BOOST_TEST(rv.end() == nullptr);
}

BOOST_AUTO_TEST_CASE(non_empty)
{
    row r = makerow("abc", 24, "def");
    row_view rv(r);

    BOOST_TEST(rv.size() == 3u);
    BOOST_TEST(rv[0] == field_view("abc"));
    BOOST_TEST(rv[1] == field_view(24));
    BOOST_TEST(rv[2] == field_view("def"));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(as_vector)
BOOST_AUTO_TEST_CASE(empty)
{
    std::vector<field> vec{field_view("abc")};
    row r;
    r.as_vector(vec);
    BOOST_TEST(vec.empty());
}

BOOST_AUTO_TEST_CASE(non_empty)
{
    std::vector<field> vec{field_view("abc")};
    row r = makerow(42u, "abc");
    r.as_vector(vec);
    BOOST_TEST(vec.size() == 2u);
    BOOST_TEST(vec[0].as_uint64() == 42u);
    BOOST_TEST(vec[1].as_string() == "abc");
}

BOOST_AUTO_TEST_CASE(return_value)
{
    auto vec = makerow(42u, "abc").as_vector();
    BOOST_TEST(vec.size() == 2u);
    BOOST_TEST(vec[0].as_uint64() == 42u);
    BOOST_TEST(vec[1].as_string() == "abc");
}
BOOST_AUTO_TEST_SUITE_END()

// operator== relies on row_view's operator==, so only
// a small subset of tests here
BOOST_AUTO_TEST_SUITE(operator_equals)
BOOST_AUTO_TEST_CASE(row_row)
{
    row r1 = makerow("abc", 4);
    row r2 = r1;
    row r3 = makerow(nullptr, 4);

    BOOST_TEST(r1 == r2);
    BOOST_TEST(!(r1 != r2));

    BOOST_TEST(!(r1 == r3));
    BOOST_TEST(r1 != r3);
}

BOOST_AUTO_TEST_CASE(row_rowview)
{
    row r1 = makerow("abc", 4);
    row r2 = makerow(nullptr, 4);
    auto fields = make_fv_arr("abc", 4);
    auto rv = makerowv(fields.data(), fields.size());

    BOOST_TEST(r1 == rv);
    BOOST_TEST(!(r1 != rv));
    BOOST_TEST(rv == r1);
    BOOST_TEST(!(rv != r1));

    BOOST_TEST(!(r2 == rv));
    BOOST_TEST(r2 != rv);
    BOOST_TEST(!(rv == r2));
    BOOST_TEST(rv != r2);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()  // test_row
