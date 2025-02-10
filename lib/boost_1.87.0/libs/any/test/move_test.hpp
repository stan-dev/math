//  See http://www.boost.org for most recent version, including documentation.
//
//  Copyright Antony Polukhin, 2013-2024.
//  Copyright Ruslan Arutyunyan, 2019.
//
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt).

#ifndef BOOST_ANY_TESTS_MOVE_TEST_HPP_INCLUDED
#define BOOST_ANY_TESTS_MOVE_TEST_HPP_INCLUDED

#include <cstdlib>
#include <string>
#include <utility>

#include <boost/core/lightweight_test.hpp>
#include <boost/type_index.hpp>

namespace any_tests {

template <typename Any>
struct move_tests // test definitions
{
    static void test_move_construction()
    {
        Any value0 = move_copy_conting_class();
        move_copy_conting_class::copy_count = 0;
        move_copy_conting_class::moves_count = 0;
        Any value(std::move(value0));

        BOOST_TEST(value0.empty());
        BOOST_TEST(!value.empty());
        BOOST_TEST(value.type() == boost::typeindex::type_id<move_copy_conting_class>());
        BOOST_TEST(boost::any_cast<move_copy_conting_class>(&value));
        BOOST_TEST_EQ(move_copy_conting_class::copy_count, 0u);
        BOOST_TEST_EQ(move_copy_conting_class::moves_count, 0u);
    }

    static void test_move_assignment()
    {
        Any value0 = move_copy_conting_class();
        Any value = move_copy_conting_class();
        move_copy_conting_class::copy_count = 0;
        move_copy_conting_class::moves_count = 0;
        value = std::move(value0);

        BOOST_TEST(value0.empty());
        BOOST_TEST(!value.empty());
        BOOST_TEST(value.type() == boost::typeindex::type_id<move_copy_conting_class>());
        BOOST_TEST(boost::any_cast<move_copy_conting_class>(&value));
        BOOST_TEST_EQ(move_copy_conting_class::copy_count, 0u);
        BOOST_TEST_EQ(move_copy_conting_class::moves_count, 0u);
    }

    static void test_copy_construction()
    {
        Any value0 = move_copy_conting_class();
        move_copy_conting_class::copy_count = 0;
        move_copy_conting_class::moves_count = 0;
        Any value(value0);

        BOOST_TEST(!value0.empty());
        BOOST_TEST(!value.empty());
        BOOST_TEST(value.type() == boost::typeindex::type_id<move_copy_conting_class>());
        BOOST_TEST(boost::any_cast<move_copy_conting_class>(&value));
        BOOST_TEST_EQ(move_copy_conting_class::copy_count, 1u);
        BOOST_TEST_EQ(move_copy_conting_class::moves_count, 0u);
    }

    static void test_copy_assignment()
    {
        Any value0 = move_copy_conting_class();
        Any value = move_copy_conting_class();
        move_copy_conting_class::copy_count = 0;
        move_copy_conting_class::moves_count = 0;
        value = value0;

        BOOST_TEST(!value0.empty());
        BOOST_TEST(!value.empty());
        BOOST_TEST(value.type() == boost::typeindex::type_id<move_copy_conting_class>());
        BOOST_TEST(boost::any_cast<move_copy_conting_class>(&value));
        BOOST_TEST_EQ(move_copy_conting_class::copy_count, 1u);
        BOOST_TEST_EQ(move_copy_conting_class::moves_count, 0u);
    }

    static void test_move_construction_from_value()
    {
        move_copy_conting_class value0;
        move_copy_conting_class::copy_count = 0;
        move_copy_conting_class::moves_count = 0;

        Any value(std::move(value0));

        BOOST_TEST(!value.empty());
        BOOST_TEST(value.type() == boost::typeindex::type_id<move_copy_conting_class>());
        BOOST_TEST(boost::any_cast<move_copy_conting_class>(&value));

        BOOST_TEST_EQ(move_copy_conting_class::copy_count, 0u);
        BOOST_TEST_EQ(move_copy_conting_class::moves_count, 1u);
     }

    static void test_move_assignment_from_value()
    {
        move_copy_conting_class value0;
        Any value;
        move_copy_conting_class::copy_count = 0;
        move_copy_conting_class::moves_count = 0;
        value = std::move(value0);

        BOOST_TEST(!value.empty());
        BOOST_TEST(value.type() == boost::typeindex::type_id<move_copy_conting_class>());
        BOOST_TEST(boost::any_cast<move_copy_conting_class>(&value));

        BOOST_TEST_EQ(move_copy_conting_class::copy_count, 0u);
        BOOST_TEST_EQ(move_copy_conting_class::moves_count, 1u);
    }

    static void test_copy_construction_from_value()
    {
        move_copy_conting_class value0;
        move_copy_conting_class::copy_count = 0;
        move_copy_conting_class::moves_count = 0;
        Any value(value0);

        BOOST_TEST(!value.empty());
        BOOST_TEST(value.type() == boost::typeindex::type_id<move_copy_conting_class>());
        BOOST_TEST(boost::any_cast<move_copy_conting_class>(&value));

        BOOST_TEST_EQ(move_copy_conting_class::copy_count, 1u);
        BOOST_TEST_EQ(move_copy_conting_class::moves_count, 0u);
     }

    static void test_copy_assignment_from_value()
    {
        move_copy_conting_class value0;
        Any value;
        move_copy_conting_class::copy_count = 0;
        move_copy_conting_class::moves_count = 0;
        value = value0;

        BOOST_TEST(!value.empty());
        BOOST_TEST(value.type() == boost::typeindex::type_id<move_copy_conting_class>());
        BOOST_TEST(boost::any_cast<move_copy_conting_class>(&value));

        BOOST_TEST_EQ(move_copy_conting_class::copy_count, 1u);
        BOOST_TEST_EQ(move_copy_conting_class::moves_count, 0u);
    }

    static const Any helper_method() {
        return true;
    }

    static const bool helper_method1() {
        return true;
    }

    static void test_construction_from_const_any_rv()
    {
        Any values[] = {helper_method(), helper_method1() };
        (void)values;
    }

    static void test_cast_to_rv()
    {
        move_copy_conting_class value0;
        Any value;
        value = value0;
        move_copy_conting_class::copy_count = 0;
        move_copy_conting_class::moves_count = 0;

        move_copy_conting_class value1 = boost::any_cast<move_copy_conting_class&&>(value);

        BOOST_TEST_EQ(move_copy_conting_class::copy_count, 0u);
        BOOST_TEST_EQ(move_copy_conting_class::moves_count, 1u);
        (void)value1;

        const Any cvalue = value0;
        move_copy_conting_class::copy_count = 0;
        move_copy_conting_class::moves_count = 0;

        move_copy_conting_class value2 = boost::any_cast<const move_copy_conting_class&>(cvalue);

        BOOST_TEST_EQ(move_copy_conting_class::copy_count, 1u);
        BOOST_TEST_EQ(move_copy_conting_class::moves_count, 0u);
        (void)value2;
    }

    class move_copy_conting_class {
    public:
        static unsigned int moves_count;
        static unsigned int copy_count;

        move_copy_conting_class(){}
        move_copy_conting_class(move_copy_conting_class&& /*param*/) {
            ++ moves_count;
        }

        move_copy_conting_class& operator=(move_copy_conting_class&& /*param*/) {
            ++ moves_count;
            return *this;
        }

        move_copy_conting_class(const move_copy_conting_class&) {
            ++ copy_count;
        }
        move_copy_conting_class& operator=(const move_copy_conting_class& /*param*/) {
            ++ copy_count;
            return *this;
        }
    };

    static int run_tests() {
        test_move_construction();
        test_move_assignment();
        test_copy_construction();
        test_copy_assignment();
        
        test_move_construction_from_value();
        test_move_assignment_from_value();
        test_copy_construction_from_value();
        test_copy_assignment_from_value();
        test_construction_from_const_any_rv();
        test_cast_to_rv();

        return boost::report_errors();
    }
}; // struct move_tests

template <typename Any>
unsigned int move_tests<Any>::move_copy_conting_class::moves_count = 0;

template <typename Any>
unsigned int move_tests<Any>::move_copy_conting_class::copy_count = 0;

} // namespace any_tests

#endif // BOOST_ANY_TESTS_MOVE_TEST_HPP_INCLUDED
