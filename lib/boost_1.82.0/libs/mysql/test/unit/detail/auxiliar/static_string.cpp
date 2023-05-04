//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/detail/auxiliar/static_string.hpp>

#include <boost/test/unit_test.hpp>

#include "printing.hpp"

using boost::mysql::string_view;
using boost::mysql::detail::static_string;

namespace {

struct test_static_string
{
    static constexpr std::size_t max_size_value = 32;
    using string_type = static_string<max_size_value>;

    static std::string original_midsize() { return "abc"; }
    static std::string original_maxsize() { return std::string(max_size_value, 'a'); }

    std::string midsize = original_midsize();
    std::string maxsize = original_maxsize();

    void wipe_midsize() { midsize = "fff"; }
    void wipe_maxsize() { maxsize = std::string(max_size_value, 'f'); }
};

// Default ctor.
BOOST_FIXTURE_TEST_CASE(default_constructor, test_static_string)
{
    string_type v;
    BOOST_TEST(v.value() == "");
}

// Init ctor.
BOOST_FIXTURE_TEST_CASE(initializing_constructor_empty_arg, test_static_string)
{
    string_type v("");
    BOOST_TEST(v.value() == "");
}

BOOST_FIXTURE_TEST_CASE(initializing_constructor_mid_size_arg, test_static_string)
{
    string_type v(midsize);
    wipe_midsize();
    BOOST_TEST(v.value() == original_midsize());
}

BOOST_FIXTURE_TEST_CASE(initializing_constructor_max_size_arg, test_static_string)
{
    string_type v(maxsize);
    wipe_maxsize();
    BOOST_TEST(v.value() == original_maxsize());
}

// Copy ctor.
BOOST_FIXTURE_TEST_CASE(copy_constructor_empty_arg, test_static_string)
{
    string_type v(string_type{});  // {} prevent deambiguation as function declaration
    BOOST_TEST(v.value() == "");
}

BOOST_FIXTURE_TEST_CASE(copy_constructor_mid_size_arg, test_static_string)
{
    string_type v(string_type{midsize});
    wipe_midsize();
    BOOST_TEST(v.value() == original_midsize());
}

BOOST_FIXTURE_TEST_CASE(copy_constructor_max_size_arg, test_static_string)
{
    string_type v(string_type{maxsize});
    wipe_maxsize();
    BOOST_TEST(v.value() == original_maxsize());
}

// Copy assignment
BOOST_FIXTURE_TEST_CASE(copy_assignment_empty_source, test_static_string)
{
    string_type v(maxsize);
    v = string_type();
    BOOST_TEST(v.value() == "");
}

BOOST_FIXTURE_TEST_CASE(copy_assignment_mid_size_source, test_static_string)
{
    string_type v(maxsize);
    v = string_type(midsize);
    wipe_midsize();
    BOOST_TEST(v.value() == original_midsize());
}

BOOST_FIXTURE_TEST_CASE(copy_assignment_max_size_source, test_static_string)
{
    string_type v(midsize);
    v = string_type(maxsize);
    wipe_midsize();
    wipe_maxsize();
    BOOST_TEST(v.value() == original_maxsize());
}

// operator==
BOOST_FIXTURE_TEST_CASE(operator_equals_both_empty, test_static_string)
{
    BOOST_TEST(string_type() == string_type());
}

BOOST_FIXTURE_TEST_CASE(operator_equals_both_empty_after_clear, test_static_string)
{
    string_type s1("abc");
    string_type s2("def");
    s1.clear();
    s2.clear();
    BOOST_TEST(s1 == s2);
}

BOOST_FIXTURE_TEST_CASE(operator_equals_one_empty_one_not, test_static_string)
{
    BOOST_TEST(!(string_type() == string_type(midsize)));
    BOOST_TEST(!(string_type(midsize) == string_type()));
    BOOST_TEST(!(string_type() == string_type(maxsize)));
    BOOST_TEST(!(string_type(maxsize) == string_type()));
}

BOOST_FIXTURE_TEST_CASE(operator_equals_same_beginning_different_size, test_static_string)
{
    string_type s1("abcd");
    string_type s2("abcde");
    BOOST_TEST(!(s1 == s2));
    BOOST_TEST(!(s2 == s1));
}

BOOST_FIXTURE_TEST_CASE(operator_equals_same_size_different_contents, test_static_string)
{
    string_type s1("abcd");
    string_type s2("dcba");
    BOOST_TEST(!(s1 == s2));
    BOOST_TEST(!(s2 == s1));
}

BOOST_FIXTURE_TEST_CASE(operator_equals_same_contents, test_static_string)
{
    BOOST_TEST(string_type(midsize) == string_type(midsize));
    BOOST_TEST(string_type(maxsize) == string_type(maxsize));
}

// operator !=
BOOST_FIXTURE_TEST_CASE(operator_not_equals_equals, test_static_string)
{
    BOOST_TEST(!(string_type() != string_type()));
    BOOST_TEST(!(string_type(midsize) != string_type(midsize)));
    BOOST_TEST(!(string_type(maxsize) != string_type(maxsize)));
}

BOOST_FIXTURE_TEST_CASE(operator_not_equals_not_equals, test_static_string)
{
    BOOST_TEST(string_type() != string_type(midsize));
    BOOST_TEST(string_type("abc") != string_type("cba"));
    BOOST_TEST(string_type(midsize) != string_type(maxsize));
}

// clear
BOOST_FIXTURE_TEST_CASE(clear_empty, test_static_string)
{
    string_type v;
    v.clear();
    BOOST_TEST(v.value() == "");
}

BOOST_FIXTURE_TEST_CASE(clear_not_empty, test_static_string)
{
    string_type v(maxsize);
    v.clear();
    BOOST_TEST(v.value() == "");
}

// append
BOOST_FIXTURE_TEST_CASE(append_from_empty_to_empty, test_static_string)
{
    string_type v;
    v.append(midsize.data(), 0);
    wipe_midsize();
    BOOST_TEST(v.value() == "");
}

BOOST_FIXTURE_TEST_CASE(append_from_empty_to_midsize, test_static_string)
{
    string_type v;
    v.append(midsize.data(), midsize.size());
    wipe_midsize();
    BOOST_TEST(v.value() == original_midsize());
}

BOOST_FIXTURE_TEST_CASE(append_from_empty_to_maxsize, test_static_string)
{
    string_type v;
    v.append(maxsize.data(), maxsize.size());
    wipe_maxsize();
    BOOST_TEST(v.value() == original_maxsize());
}

BOOST_FIXTURE_TEST_CASE(append_from_midsize_to_midsize, test_static_string)
{
    string_type v("222");
    v.append(midsize.data(), midsize.size());
    wipe_midsize();
    BOOST_TEST(v.value() == "222" + original_midsize());
}

BOOST_FIXTURE_TEST_CASE(append_from_midsize_to_maxsize, test_static_string)
{
    string_type v(midsize);
    std::string newbuff(max_size_value - midsize.size(), '1');
    v.append(newbuff.data(), newbuff.size());
    wipe_midsize();
    BOOST_TEST(v.value() == original_midsize() + newbuff);
}

}  // namespace
