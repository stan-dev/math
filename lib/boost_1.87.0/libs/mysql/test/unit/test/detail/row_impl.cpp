//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/blob.hpp>
#include <boost/mysql/blob_view.hpp>
#include <boost/mysql/date.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/row_view.hpp>

#include <boost/mysql/detail/row_impl.hpp>

#include <boost/test/unit_test.hpp>

#include <algorithm>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/create_basic.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;
using boost::mysql::detail::row_impl;

BOOST_AUTO_TEST_SUITE(test_row_impl)

template <class... T>
row_impl makerowimpl(T&&... args)
{
    auto fields = make_fv_arr(std::forward<T>(args)...);
    return row_impl(fields.data(), fields.size());
}

template <class... T>
void add_fields(row_impl& r, T&&... args)
{
    auto fields = make_fv_arr(std::forward<T>(args)...);
    std::copy(fields.begin(), fields.end(), r.add_fields(fields.size()).data());
}

// An array and vector with all scalar types
std::vector<field_view> make_scalar_vector()
{
    return make_fv_vector(
        42,
        42u,
        4.2f,
        4.2,
        date(2020, 1, 2),
        datetime(2020, 10, 10),
        maket(10, 1, 1),
        nullptr
    );
}

void hard_clear(std::vector<field_view>& res)
{
    for (auto& f : res)
        f = field_view();
}

// Check that references to certain row_impl's contents survive an operation
struct reference_checker
{
    const field_view* ptr;

    reference_checker(const row_impl& r) : ptr(r.fields().data()) {}

    void check(const row_impl& new_row) { BOOST_TEST(ptr == new_row.fields().data()); }
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

    reference_checker_strs(const row_impl& r, std::size_t string_index, std::size_t blob_index)
        : reference_checker(r),
          string_index(string_index),
          blob_index(blob_index),
          string_ptr(r.fields().at(string_index).as_string().data()),
          blob_ptr(r.fields().at(blob_index).as_blob().data())
    {
    }

    void check(const row_impl& new_row)
    {
        reference_checker::check(new_row);
        assert_ptrs_equal(string_ptr, new_row.fields().at(string_index).as_string().data());
        BOOST_TEST(blob_ptr == new_row.fields().at(blob_index).as_blob().data());
    }
};

BOOST_AUTO_TEST_CASE(default_ctor)
{
    row_impl r;
    BOOST_TEST(r.fields().empty());
}

BOOST_AUTO_TEST_SUITE(ctor_from_span)
BOOST_AUTO_TEST_CASE(empty)
{
    std::vector<field_view> fields;
    row_impl r(fields.data(), fields.size());
    BOOST_TEST(r.fields().empty());
}

BOOST_AUTO_TEST_CASE(non_strings)
{
    auto fields = make_scalar_vector();
    row_impl r(fields.data(), fields.size());

    // Fields still valid even when the original source of the view changed
    hard_clear(fields);
    BOOST_TEST(r.fields() == make_scalar_vector());
}

BOOST_AUTO_TEST_CASE(strings_blobs)
{
    std::string s1("test"), s2("othertest");
    blob b{0x00, 0xab, 0xf5};
    auto fields = make_fv_arr(s1, s2, 50, b);
    row_impl r(fields.data(), fields.size());

    // Fields still valid even when the original strings changed
    s1 = "jkiop";
    s2 = "abcdef";
    b = {0xff, 0xa4, 0x02};
    BOOST_TEST(r.fields().size() == 4u);
    BOOST_TEST(r.fields()[0] == field_view("test"));
    BOOST_TEST(r.fields()[1] == field_view("othertest"));
    BOOST_TEST(r.fields()[2] == field_view(50));
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(r.fields()[3].as_blob(), blob({0x00, 0xab, 0xf5}));
}

BOOST_AUTO_TEST_CASE(empty_strings_blobs)
{
    std::string s;
    blob b;
    auto fields = make_fv_arr(s, 50, b);
    row_impl r(fields.data(), fields.size());

    // Fields still valid even when the original strings changed
    s = "other";
    b = {0xff, 0xa4, 0x02};
    BOOST_TEST(r.fields().size() == 3u);
    BOOST_TEST(r.fields()[0] == field_view(""));
    BOOST_TEST(r.fields()[1] == field_view(50));
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(r.fields()[2].as_blob(), blob());
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(copy_ctor)
BOOST_AUTO_TEST_CASE(empty)
{
    row_impl r1;
    row_impl r2(r1);
    r1 = makerowimpl(42, "test");  // r2 should be independent of r1

    BOOST_TEST(r2.fields().empty());
}

BOOST_AUTO_TEST_CASE(non_strings)
{
    auto fields = make_scalar_vector();
    row_impl r1(fields.data(), fields.size());
    row_impl r2(r1);
    r1 = makerowimpl(42, "test");  // r2 should be independent of r1

    BOOST_TEST(r2.fields() == make_scalar_vector());
}

BOOST_AUTO_TEST_CASE(strings_blobs)
{
    row_impl r1 = makerowimpl("", 42, "test", makebv("\0\3\2"));
    row_impl r2(r1);
    r1 = makerowimpl("another_string", 4.2f, "");  // r2 should be independent of r1

    BOOST_TEST(r2.fields().size() == 4u);
    BOOST_TEST(r2.fields()[0] == field_view(""));
    BOOST_TEST(r2.fields()[1] == field_view(42));
    BOOST_TEST(r2.fields()[2] == field_view("test"));
    BOOST_TEST(r2.fields()[3] == field_view(makebv("\0\3\2")));
}

BOOST_AUTO_TEST_CASE(empty_strings_blobs)
{
    row_impl r1 = makerowimpl("", 42, blob_view());
    row_impl r2(r1);
    r1 = makerowimpl("another_string", 4.2f, "");  // r2 should be independent of r1

    BOOST_TEST(r2.fields().size() == 3u);
    BOOST_TEST(r2.fields()[0] == field_view(""));
    BOOST_TEST(r2.fields()[1] == field_view(42));
    BOOST_TEST(r2.fields()[2] == field_view(blob_view()));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(move_ctor)
BOOST_AUTO_TEST_CASE(empty)
{
    row_impl r1;

    // References should remain valid
    reference_checker refcheck(r1);

    row_impl r2(std::move(r1));
    r1 = makerowimpl(42, "test");  // r2 should be independent of r1

    BOOST_TEST(r2.fields().empty());
    refcheck.check(r2);
}

BOOST_AUTO_TEST_CASE(non_strings)
{
    auto fields = make_scalar_vector();
    row_impl r1(fields.data(), fields.size());

    // References, pointers, etc. should remain valid
    reference_checker refcheck(r1);

    row_impl r2(std::move(r1));
    r1 = makerowimpl(42, "test");  // r2 should be independent of r1

    BOOST_TEST(r2.fields() == make_scalar_vector());
    refcheck.check(r2);
}

BOOST_AUTO_TEST_CASE(strings_blobs)
{
    row_impl r1 = makerowimpl("", 42, "test", makebv("\0\5\xff"));

    // References, pointers, etc should remain valid
    reference_checker_strs refcheck(r1, 2, 3);

    // Move
    row_impl r2(std::move(r1));
    r1 = makerowimpl("another_string", 4.2f, "", makebv("\1\5\xab"));  // r2 should be independent of r1

    BOOST_TEST(r2.fields().size() == 4u);
    BOOST_TEST(r2.fields()[0] == field_view(""));
    BOOST_TEST(r2.fields()[1] == field_view(42));
    BOOST_TEST(r2.fields()[2] == field_view("test"));
    BOOST_TEST(r2.fields()[3] == field_view(makebv("\0\5\xff")));

    // References, pointers, etc still valid
    refcheck.check(r2);
}

BOOST_AUTO_TEST_CASE(empty_strings_blobs)
{
    row_impl r1 = makerowimpl("", 42, blob_view());

    // Move
    row_impl r2(std::move(r1));
    r1 = makerowimpl("another_string", 4.2f, "", makebv("\1\5\xab"));  // r2 should be independent of r1

    BOOST_TEST(r2.fields() == make_fv_vector("", 42, blob_view()));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(copy_assignment)
BOOST_AUTO_TEST_CASE(empty)
{
    row_impl r1 = makerowimpl(42, "abcdef");
    row_impl r2;
    r1 = r2;
    r2 = makerowimpl(90, nullptr);  // r1 is independent of r2
    BOOST_TEST(r1.fields().empty());
}

BOOST_AUTO_TEST_CASE(non_strings)
{
    auto fields = make_scalar_vector();
    row_impl r1 = makerowimpl(42, "abcdef");
    row_impl r2(fields.data(), fields.size());
    r1 = r2;
    r2 = makerowimpl("abc", 80, nullptr);  // r1 is independent of r2

    BOOST_TEST(r1.fields() == make_scalar_vector());
}

BOOST_AUTO_TEST_CASE(strings_blobs)
{
    row_impl r1 = makerowimpl(42, "abcdef", makebv("\0\1\2"));
    row_impl r2 = makerowimpl("a_very_long_string", nullptr, "", makebv("\3\4\5"));
    r1 = r2;
    r2 = makerowimpl("another_string", 90, "yet_another");  // r1 is independent of r2

    BOOST_TEST(r1.fields() == make_fv_vector("a_very_long_string", nullptr, "", makebv("\3\4\5")));
}

BOOST_AUTO_TEST_CASE(empty_strings_blobs)
{
    row_impl r1 = makerowimpl(42, "abcdef", makebv("\0\1\2"));
    row_impl r2 = makerowimpl(nullptr, "", blob_view());
    r1 = r2;
    r2 = makerowimpl("another_string", 90, "yet_another");  // r1 is independent of r2

    BOOST_TEST(r1.fields() == make_fv_vector(nullptr, "", blob_view()));
}

BOOST_AUTO_TEST_CASE(strings_blobs_empty_to)
{
    row_impl r1;
    row_impl r2 = makerowimpl("abc", nullptr, "bcd", makebv("\1\2\3"));
    r1 = r2;

    BOOST_TEST(r1.fields() == make_fv_vector("abc", nullptr, "bcd", makebv("\1\2\3")));
}

BOOST_AUTO_TEST_CASE(self_assignment_empty)
{
    row_impl r;
    const row_impl& ref = r;
    r = ref;

    BOOST_TEST(r.fields().empty());
}

BOOST_AUTO_TEST_CASE(self_assignment_non_empty)
{
    row_impl r = makerowimpl("abc", 50u, "fgh");
    const row_impl& ref = r;
    r = ref;

    BOOST_TEST(r.fields() == make_fv_vector("abc", 50u, "fgh"));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(move_assignment)
BOOST_AUTO_TEST_CASE(empty)
{
    row_impl r1 = makerowimpl(42, "abcdef");
    row_impl r2;

    reference_checker refcheck(r2);

    r1 = std::move(r2);
    r2 = makerowimpl(90, nullptr);  // r1 is independent of r2
    BOOST_TEST(r1.fields().empty());
    refcheck.check(r1);
}

BOOST_AUTO_TEST_CASE(non_strings)
{
    auto fields = make_scalar_vector();
    row_impl r1 = makerowimpl(42, "abcdef");
    row_impl r2(fields.data(), fields.size());

    // References should remain valid
    reference_checker refcheck(r2);

    r1 = std::move(r2);
    r2 = makerowimpl("abc", 80, nullptr);  // r1 is independent of r2

    BOOST_TEST(r1.fields() == make_scalar_vector());
    refcheck.check(r1);
}

BOOST_AUTO_TEST_CASE(strings_blobs)
{
    row_impl r1 = makerowimpl(42, "abcdef", makebv("\0\4\1"));
    row_impl r2 = makerowimpl("a_very_long_string", nullptr, "", makebv("\7\1\2"));

    // References, pointers, etc should remain valid
    reference_checker_strs refcheck(r2, 0, 3);

    // Move
    r1 = std::move(r2);
    r2 = makerowimpl("another_string", 90, "yet_another", makebv("\0\0"));  // r1 is independent of r2

    BOOST_TEST(r1.fields() == make_fv_vector("a_very_long_string", nullptr, "", makebv("\7\1\2")));
    refcheck.check(r1);
}

BOOST_AUTO_TEST_CASE(empty_strings_blobs)
{
    row_impl r1 = makerowimpl(42, "abcdef", makebv("\0\4\1"));
    row_impl r2 = makerowimpl("", blob_view());

    // References, pointers, etc should remain valid
    reference_checker_strs refcheck(r2, 0, 1);

    // Move
    r1 = std::move(r2);
    r2 = makerowimpl("another_string", 90);  // r1 is independent of r2

    BOOST_TEST(r1.fields() == make_fv_vector("", blob_view()));
    refcheck.check(r1);
}

BOOST_AUTO_TEST_CASE(strings_blobs_empty_to)
{
    row_impl r1;
    row_impl r2 = makerowimpl("abc", nullptr, "bcd", makebv("\0\2\5"));

    // References, pointers, etc should remain valid
    reference_checker_strs refcheck(r2, 2, 3);

    r1 = std::move(r2);

    BOOST_TEST(r1.fields() == make_fv_vector("abc", nullptr, "bcd", makebv("\0\2\5")));
    refcheck.check(r1);
}

BOOST_AUTO_TEST_CASE(self_assignment_empty)
{
    row_impl r;
    row_impl&& ref = std::move(r);
    r = std::move(ref);

    // r is in a valid but unspecified state; can be assigned to
    r = makerowimpl("abcdef");
    BOOST_TEST(r.fields() == make_fv_vector("abcdef"));
}

BOOST_AUTO_TEST_CASE(self_assignment_non_empty)
{
    row_impl r = makerowimpl("abc", 50u, "fgh", makebv("\0\4"));
    row_impl&& ref = std::move(r);
    r = std::move(ref);  // this should leave r in a valid but unspecified state

    // r is in a valid but unspecified state; can be assigned to
    r = makerowimpl("abcdef");
    BOOST_TEST(r.fields() == make_fv_vector("abcdef"));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(assignment_from_span)
BOOST_AUTO_TEST_CASE(empty)
{
    row_impl r = makerowimpl(42, "abcdef", makebv("\0\xae"));
    r.assign(nullptr, 0);
    BOOST_TEST(r.fields().empty());
}

BOOST_AUTO_TEST_CASE(empty_non_null)
{
    field_view f;
    row_impl r = makerowimpl(42, "abcdef", makebv("\0\xae"));
    r.assign(&f, 0);
    BOOST_TEST(r.fields().empty());
}

BOOST_AUTO_TEST_CASE(non_strings)
{
    row_impl r = makerowimpl(42, "abcdef");
    auto fields = make_scalar_vector();
    r.assign(fields.data(), fields.size());
    hard_clear(fields);  // r should be independent of the original fields

    BOOST_TEST(r.fields() == make_scalar_vector());
}

BOOST_AUTO_TEST_CASE(strings_blobs)
{
    std::string s1("a_very_long_string"), s2("abc");
    blob b{0x00, 0xfa};
    row_impl r = makerowimpl(42, "haksj", makebv("\0\1"));
    auto fields = make_fv_arr(s1, nullptr, s2, b);

    r.assign(fields.data(), fields.size());
    fields = make_fv_arr("abc", 42u, 9, nullptr);  // r should be independent of the original fields
    s1 = "another_string";                         // r should be independent of the original strings
    s2 = "yet_another";
    b = {0xac, 0x32, 0x21, 0x50};

    BOOST_TEST(r.fields() == make_fv_vector("a_very_long_string", nullptr, "abc", makebv("\0\xfa")));
}

BOOST_AUTO_TEST_CASE(empty_strings_blobs)
{
    std::string s("");
    blob b;
    row_impl r = makerowimpl(42, "haksj", makebv("\0\1"));
    auto fields = make_fv_arr(s, b);

    r.assign(fields.data(), fields.size());
    fields = make_fv_arr(0, 0);  // r should be independent of the original fields
    s = "another_string";        // r should be independent of the original strings
    b = {0xac, 0x32, 0x21, 0x50};

    BOOST_TEST(r.fields() == make_fv_vector("", blob_view()));
}

BOOST_AUTO_TEST_CASE(strings_blobs_empty_to)
{
    row_impl r;
    auto fields = make_fv_arr("abc", nullptr, "bcd", makebv("\0\3"));
    r.assign(fields.data(), fields.size());

    BOOST_TEST(r.fields() == make_fv_vector("abc", nullptr, "bcd", makebv("\0\3")));
}

BOOST_AUTO_TEST_CASE(self_assignment)
{
    row_impl r = makerowimpl("abcdef", 42, "plk", makebv("\0\1"));
    r.assign(r.fields().data(), r.fields().size());

    BOOST_TEST(r.fields() == make_fv_vector("abcdef", 42, "plk", makebv("\0\1")));
}

BOOST_AUTO_TEST_CASE(self_assignment_empty)
{
    row_impl r;
    r.assign(r.fields().data(), r.fields().size());
    BOOST_TEST(r.fields().empty());
}

BOOST_AUTO_TEST_CASE(self_assignment_cleared)
{
    row_impl r = makerowimpl(42, "abc");
    r.clear();
    r.assign(r.fields().data(), r.fields().size());
    BOOST_TEST(r.fields().empty());
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(add_fields_)
BOOST_AUTO_TEST_CASE(empty_collection)
{
    row_impl r;
    auto storage = r.add_fields(2);
    BOOST_TEST(r.fields().size() == 2u);
    BOOST_TEST(storage.data() == r.fields().data());
    BOOST_TEST(storage.size() == 2u);
}

BOOST_AUTO_TEST_CASE(non_empty_collection)
{
    row_impl r = makerowimpl(nullptr, nullptr);
    auto storage = r.add_fields(3);
    BOOST_TEST(r.fields().size() == 5u);
    BOOST_TEST(storage.data() == r.fields().data() + 2);
    BOOST_TEST(storage.size() == 3u);
}

BOOST_AUTO_TEST_CASE(zero_fields)
{
    row_impl r = makerowimpl(nullptr, nullptr);
    auto storage = r.add_fields(0);
    BOOST_TEST(r.fields().size() == 2u);
    BOOST_TEST(storage.data() == r.fields().data() + 2);
    BOOST_TEST(storage.size() == 0u);
}

BOOST_AUTO_TEST_CASE(empty_collection_zero_fields)
{
    row_impl r;
    auto storage = r.add_fields(0);
    BOOST_TEST(r.fields().size() == 0u);
    BOOST_TEST(storage.data() == r.fields().data());
    BOOST_TEST(storage.size() == 0u);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(copy_strings_as_offsets)
BOOST_AUTO_TEST_CASE(scalars)
{
    row_impl r;
    add_fields(r, nullptr, 42, 10.0f, date(2020, 10, 1));
    r.copy_strings_as_offsets(0, 4);
    r.offsets_to_string_views();
    BOOST_TEST(r.fields() == make_fv_vector(nullptr, 42, 10.f, date(2020, 10, 1)));
}

BOOST_AUTO_TEST_CASE(strings_blobs)
{
    row_impl r;
    std::string s = "abc";
    blob b{0x01, 0x02, 0x03};
    add_fields(r, nullptr, s, 10.f, b);
    r.copy_strings_as_offsets(1, 3);
    s = "ghi";
    b = {0xff, 0xff, 0xff};
    r.offsets_to_string_views();
    BOOST_TEST(r.fields() == make_fv_vector(nullptr, "abc", 10.f, makebv("\1\2\3")));
}

BOOST_AUTO_TEST_CASE(empty_strings_blobs)
{
    row_impl r;
    std::string s = "";
    blob b{};
    add_fields(r, nullptr, s, 10.f, b);
    r.copy_strings_as_offsets(1, 3);
    s = "ghi";
    b = {0xff, 0xff, 0xff};
    r.offsets_to_string_views();
    BOOST_TEST(r.fields() == make_fv_vector(nullptr, "", 10.f, makebv("")));
}

BOOST_AUTO_TEST_CASE(buffer_relocation)
{
    row_impl r;
    std::string s = "abc";
    add_fields(r, nullptr, s);
    r.copy_strings_as_offsets(0, 2);
    s = "ghi";

    blob b{0x01, 0x02, 0x03};
    add_fields(r, 10.f, b);
    r.copy_strings_as_offsets(2, 2);

    s = "";
    b = {};
    add_fields(r, s, b);
    r.copy_strings_as_offsets(4, 2);
    b = {0x01, 0x02};

    s = "this is a long string";
    add_fields(r, s);
    r.copy_strings_as_offsets(6, 1);
    s = "another long string";

    r.offsets_to_string_views();
    BOOST_TEST(
        r.fields() ==
        make_fv_vector(nullptr, "abc", 10.f, makebv("\1\2\3"), "", makebv(""), "this is a long string")
    );
}

BOOST_AUTO_TEST_CASE(empty_range)
{
    std::string s = "abc";
    row_impl r = makerowimpl(nullptr, 42);
    r.copy_strings_as_offsets(0, 0);
    r.offsets_to_string_views();
    BOOST_TEST(r.fields() == make_fv_vector(nullptr, 42));
}

BOOST_AUTO_TEST_CASE(empty_collection)
{
    row_impl r;
    r.copy_strings_as_offsets(0, 0);
    r.offsets_to_string_views();
    BOOST_TEST(r.fields().empty());
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
