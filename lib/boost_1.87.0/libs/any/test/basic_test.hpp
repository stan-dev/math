// Copyright Kevlin Henney, 2000, 2001. All rights reserved.
// Copyright Antony Polukhin, 2013-2019.
// Copyright Ruslan Arutyunyan, 2019-2021.
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_ANY_TEST_BASIC_TEST_HPP
#define BOOST_ANY_TEST_BASIC_TEST_HPP

// what:  unit tests for variant type boost::any
// who:   contributed by Kevlin Henney
// when:  July 2001, 2013, 2014
// where: tested with BCC 5.5, MSVC 6.0, and g++ 2.95

#include <cstdlib>
#include <string>
#include <vector>
#include <utility>
#include <type_traits>

#include <boost/core/lightweight_test.hpp>

namespace any_tests {

struct huge_structure {
    char take_place[1024];
    std::string text;
};


template <typename Any>
struct basic_tests  // test definitions
{
    struct copy_counter
    {
    public:
        copy_counter() {}
        copy_counter(const copy_counter&) { ++count; }
        copy_counter& operator=(const copy_counter&) { ++count; return *this; }
        static int get_count() { return count; }

    private:
        static int count;
    };

    static void test_default_ctor()
    {
        const Any value;

        BOOST_TEST(value.empty());
        BOOST_TEST_EQ(boost::any_cast<int>(&value), nullptr);
        BOOST_TEST(value.type() == boost::typeindex::type_id<void>());
    }

    static void test_converting_ctor()
    {
        std::string text = "test message";
        Any value = text;

        BOOST_TEST(!value.empty());
        BOOST_TEST(value.type() == boost::typeindex::type_id<std::string>());
        BOOST_TEST_EQ(boost::any_cast<int>(&value), nullptr);
        BOOST_TEST_NE(boost::any_cast<std::string>(&value), nullptr);
        BOOST_TEST_EQ(boost::any_cast<std::string>(value), text);
        BOOST_TEST_NE(boost::any_cast<std::string>(&value), &text);
    }

    static void test_copy_ctor()
    {
        std::string text = "test message";
        Any original = text, copy = original;

        BOOST_TEST(!copy.empty());
        BOOST_TEST(boost::typeindex::type_index(original.type()) == copy.type());
        BOOST_TEST_EQ(
            boost::any_cast<std::string>(original),
            boost::any_cast<std::string>(copy)
        );
        BOOST_TEST_EQ(text, boost::any_cast<std::string>(copy));
        BOOST_TEST_NE(
            boost::any_cast<std::string>(&original),
            boost::any_cast<std::string>(&copy)
        );
    }

    static void test_copy_assign()
    {
        std::string text = "test message";
        Any original = text, copy;
        Any * assign_result = &(copy = original);

        BOOST_TEST(!copy.empty());
        BOOST_TEST(boost::typeindex::type_index(original.type()) == copy.type());
        BOOST_TEST_EQ(
            boost::any_cast<std::string>(original),
            boost::any_cast<std::string>(copy)
        );
        BOOST_TEST_EQ(text, boost::any_cast<std::string>(copy));
        BOOST_TEST_NE(
            boost::any_cast<std::string>(&original),
            boost::any_cast<std::string>(&copy)
        );
        BOOST_TEST_EQ(assign_result, &copy);
    }

    static void test_converting_assign()
    {
        std::string text = "test message";
        Any value;
        Any * assign_result = &(value = text);

        BOOST_TEST(!value.empty());
        BOOST_TEST(value.type() == boost::typeindex::type_id<std::string>());
        BOOST_TEST_EQ(boost::any_cast<int>(&value), nullptr);
        BOOST_TEST_NE(boost::any_cast<std::string>(&value), nullptr);
        BOOST_TEST_EQ(boost::any_cast<std::string>(value), text);
        BOOST_TEST_NE(boost::any_cast<std::string>(&value), &text);
        BOOST_TEST_EQ(assign_result, &value);
    }

    static void test_bad_cast()
    {
        std::string text = "test message";
        Any value = text;

        BOOST_TEST_THROWS(
            boost::any_cast<const char *>(value),
            boost::bad_any_cast
        );
    }

    static void test_swap()
    {
        huge_structure stored;
        stored.text = "test message";

        Any original = stored;
        Any swapped;
        huge_structure * original_ptr = boost::any_cast<huge_structure>(&original);
        Any * swap_result = &original.swap(swapped);

        BOOST_TEST(original.empty());
        BOOST_TEST(!swapped.empty());
        BOOST_TEST(swapped.type() == boost::typeindex::type_id<huge_structure>());
        BOOST_TEST_EQ(stored.text, boost::any_cast<huge_structure>(swapped).text);
        BOOST_TEST_NE(original_ptr, nullptr);
        BOOST_TEST_EQ(
            original_ptr,
            boost::any_cast<huge_structure>(&swapped)
        );
        BOOST_TEST_EQ(swap_result, &original);

        swap(swapped, swapped);
        BOOST_TEST(!swapped.empty());
        BOOST_TEST(
            swapped.type() == boost::typeindex::type_id<huge_structure>()
        );
        BOOST_TEST_EQ(
            stored.text, boost::any_cast<huge_structure>(swapped).text
        );

        Any copy1 = copy_counter();
        Any copy2 = copy_counter();
        int count = copy_counter::get_count();
        swap(copy1, copy2);
        BOOST_TEST_EQ(count, copy_counter::get_count());

        Any any_char = '1';
        swap(any_char, swapped);
        BOOST_TEST_EQ(
            stored.text, boost::any_cast<huge_structure>(any_char).text
        );
        BOOST_TEST(
            swapped.type() == boost::typeindex::type_id<char>()
        );
        BOOST_TEST_EQ('1', boost::any_cast<char>(swapped));
    }

    static void test_null_copying()
    {
        const Any null;
        Any copied = null, assigned;
        assigned = null;

        BOOST_TEST(null.empty());
        BOOST_TEST(copied.empty());
        BOOST_TEST(assigned.empty());
    }

    static void test_cast_to_reference()
    {
        Any a(137);
        const Any b(a);

        int &                ra    = boost::any_cast<int &>(a);
        int const &          ra_c  = boost::any_cast<int const &>(a);
        int volatile &       ra_v  = boost::any_cast<int volatile &>(a);
        int const volatile & ra_cv = boost::any_cast<int const volatile&>(a);

        BOOST_TEST(&ra == &ra_c && &ra == &ra_v && &ra == &ra_cv);

        int const &          rb_c  = boost::any_cast<int const &>(b);
        int const volatile & rb_cv = boost::any_cast<int const volatile &>(b);

        BOOST_TEST(&rb_c == &rb_cv);
        BOOST_TEST(&ra != &rb_c);

        ++ra;
        int incremented = boost::any_cast<int>(a);
        BOOST_TEST(incremented == 138);

        BOOST_TEST_THROWS(
            boost::any_cast<char &>(a),
            boost::bad_any_cast
        );

        BOOST_TEST_THROWS(
            boost::any_cast<const char &>(b),
            boost::bad_any_cast
        );
    }

    static void test_bad_any_cast()
    {
        BOOST_TEST((
            std::is_base_of<std::exception, boost::bad_any_cast>::value
        ));

        BOOST_TEST(
            std::string(boost::bad_any_cast().what()).find("any") != std::string::npos
        );
    }

    static void test_with_array()
    {
        Any value1("Char array");
        Any value2;
        value2 = "Char array";

        BOOST_TEST(!value1.empty());
        BOOST_TEST(!value2.empty());

        BOOST_TEST(value1.type() == boost::typeindex::type_id<const char*>());
        BOOST_TEST(value2.type() == boost::typeindex::type_id<const char*>());

        BOOST_TEST(boost::any_cast<const char*>(&value1));
        BOOST_TEST(boost::any_cast<const char*>(&value2));
    }

    static void test_clear()
    {
        std::string text = "test message";
        Any value = text;

        BOOST_TEST(!value.empty());

        value.clear();
        BOOST_TEST(value.empty());

        value.clear();
        BOOST_TEST(value.empty());

        value = text;
        BOOST_TEST(!value.empty());

        value.clear();
        BOOST_TEST(value.empty());
    }

    // Following tests cover the case from #9462
    // https://svn.boost.org/trac/boost/ticket/9462
    static Any makeVec()
    {
        return std::vector<int>(100 /*size*/, 7 /*value*/);
    }

    static void test_vectors()
    {
        const std::vector<int>& vec = boost::any_cast<std::vector<int> >(makeVec());
        BOOST_TEST_EQ(vec.size(), 100u);
        BOOST_TEST_EQ(vec.back(), 7);
        BOOST_TEST_EQ(vec.front(), 7);

        std::vector<int> vec1 = boost::any_cast<std::vector<int> >(makeVec());
        BOOST_TEST_EQ(vec1.size(), 100u);
        BOOST_TEST_EQ(vec1.back(), 7);
        BOOST_TEST_EQ(vec1.front(), 7);
    }

    template<typename T>
    class class_with_address_op
    {
    public:
        class_with_address_op(const T* p)
            : ptr(p)
        {}

        const T** operator &()
        {
            return &ptr;
        }

        const T* get() const
        {
            return ptr;
        }

    private:
        const T* ptr;
    };

    static void test_addressof()
    {
        int val = 10;
        const int* ptr = &val;
        class_with_address_op<int> obj(ptr);
        Any test_val(obj);

        class_with_address_op<int> returned_obj = boost::any_cast<class_with_address_op<int> >(test_val);
        BOOST_TEST_EQ(&val, returned_obj.get());

        BOOST_TEST(!!boost::any_cast<class_with_address_op<int> >(&test_val));
        BOOST_TEST_EQ(boost::unsafe_any_cast<class_with_address_op<int> >(&test_val)->get(), ptr);
    }

    static void test_multiple_assign()
    {
        Any test_val = 10;
        BOOST_TEST(!!boost::any_cast<int>(&test_val));

        test_val = '0';
        BOOST_TEST(!!boost::any_cast<char>(&test_val));

        test_val = huge_structure();
        BOOST_TEST(!!boost::any_cast<huge_structure>(&test_val));

        test_val = '0';
        BOOST_TEST(!!boost::any_cast<char>(&test_val));

        test_val = Any(huge_structure());
        BOOST_TEST(!!boost::any_cast<huge_structure>(&test_val));
    }

    static int run_tests()
    {
        test_default_ctor();
        test_converting_ctor();
        test_copy_ctor();
        test_copy_assign();
        test_converting_assign();
        test_bad_cast();
        test_swap();
        test_null_copying();
        test_cast_to_reference();
        test_bad_any_cast();
        test_with_array();
        test_clear();
        test_vectors();
        test_addressof();
        test_multiple_assign();

        return boost::report_errors();
    }
};


template <typename Any>
int basic_tests<Any>::copy_counter::count = 0;

}

#endif
