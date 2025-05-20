// Copyright Antony Polukhin, 2013-2024.
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/any/unique_any.hpp>

#include <boost/core/lightweight_test.hpp>

void test_move_construct_unique_ptr() {
    std::unique_ptr<int> ptr(new int(42));
    auto* raw_ptr = ptr.get();

    boost::anys::unique_any a = std::move(ptr);
    BOOST_TEST(!ptr);
    BOOST_TEST(a.has_value());
    BOOST_TEST_EQ(boost::any_cast<std::unique_ptr<int>&>(a).get(), raw_ptr);
}

void test_move_construct_unique_any() {
    std::unique_ptr<int> ptr(new int(42));
    auto* raw_ptr = ptr.get();

    boost::anys::unique_any a = std::move(ptr);
    boost::anys::unique_any b = std::move(a);
    BOOST_TEST(!a.has_value());
    BOOST_TEST_EQ(boost::any_cast<std::unique_ptr<int>&>(b).get(), raw_ptr);
}

void test_move_assign_unique_ptr() {
    std::unique_ptr<int> ptr(new int(42));
    auto* raw_ptr = ptr.get();

    boost::anys::unique_any a;
    a = std::move(ptr);
    BOOST_TEST(!ptr);
    BOOST_TEST(a.has_value());
    BOOST_TEST_EQ(boost::any_cast<std::unique_ptr<int>&>(a).get(), raw_ptr);
}

void test_move_assign_unique_any() {
    std::unique_ptr<int> ptr(new int(42));
    auto* raw_ptr = ptr.get();

    boost::anys::unique_any a = std::move(ptr);
    boost::anys::unique_any b;
    b = std::move(a);
    BOOST_TEST(!a.has_value());
    BOOST_TEST_EQ(boost::any_cast<std::unique_ptr<int>&>(b).get(), raw_ptr);
}

void test_move_any_cast_implicit() {
    std::unique_ptr<int> ptr(new int(42));
    auto* raw_ptr = ptr.get();

    boost::anys::unique_any a = std::move(ptr);
    auto new_ptr = boost::any_cast<std::unique_ptr<int>>(std::move(a));
    BOOST_TEST_EQ(new_ptr.get(), raw_ptr);
}

void test_move_any_cast_explicit() {
    std::unique_ptr<int> ptr(new int(42));
    auto* raw_ptr = ptr.get();

    boost::anys::unique_any a = std::move(ptr);
    auto new_ptr = boost::any_cast<std::unique_ptr<int>&&>(std::move(a));
    BOOST_TEST_EQ(new_ptr.get(), raw_ptr);
}

int main() {
    test_move_construct_unique_ptr();
    test_move_construct_unique_any();
    test_move_assign_unique_ptr();
    test_move_assign_unique_any();
    test_move_any_cast_implicit();
    test_move_any_cast_explicit();

    return boost::report_errors();
}

