//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
// Copyright (c) 2022 Alexander Grund
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include "../src/boost/locale/shared/ios_prop.hpp"
#include "boostLocale/test/unit_test.hpp"
#include <locale>
#include <sstream>

int counter = 0;
int imbued = 0;
struct test_property {
    explicit test_property(int xx = -1) : x(xx) { counter++; }
    test_property(const test_property& other)
    {
        counter++;
        x = other.x;
    }
    test_property& operator=(const test_property&) = default;
    ~test_property() { counter--; }
    void on_imbue() { imbued++; }

    int x;
};
typedef boost::locale::impl::ios_prop<test_property> prop_type;

struct init {
    init() { prop_type::global_init(); }
};

void test_main(int /*argc*/, char** /*argv*/)
{
    // Test get() default constructs the property if it does not exist or returns the existing one
    {
        std::stringstream ss;
        TEST_EQ(prop_type::get(ss).x, -1);
        TEST_EQ(counter, 1);
        prop_type::get(ss).x = 42;
        TEST_EQ(prop_type::get(ss).x, 42);
        TEST_EQ(counter, 1);
    }
    TEST_EQ(counter, 0);
    {
        std::stringstream ss;
        prop_type::get(ss).x = 1;
        {
            std::stringstream ss2;
            TEST_EQ(counter, 1);
            ss2.copyfmt(ss);
            TEST_EQ(counter, 2);
            TEST_EQ(prop_type::get(ss).x, 1);
            TEST_EQ(prop_type::get(ss2).x, 1);
            // Check that both are distinct copies
            prop_type::get(ss2).x = 2;
            TEST_EQ(prop_type::get(ss).x, 1);
            TEST_EQ(prop_type::get(ss2).x, 2); // Only 2nd is changed
            // Imbue on 2nd calls on_imbue once
            TEST_EQ(imbued, 0);
            ss2.imbue(std::locale::classic());
            TEST_EQ(imbued, 1);
            // Copy again over existing
            ss2.copyfmt(ss);
            TEST_EQ(counter, 2);
            TEST_EQ(prop_type::get(ss).x, 1);
            TEST_EQ(prop_type::get(ss2).x, 1); // Copied
            prop_type::get(ss2).x = 2;
            TEST_EQ(prop_type::get(ss).x, 1);
            TEST_EQ(prop_type::get(ss2).x, 2); // Only this changed
        }
        // Copy from unset removes the property
        {
            std::stringstream ss2;
            TEST_EQ(prop_type::get(ss).x, 1);
            TEST_EQ(counter, 1);
            ss.copyfmt(ss2);
            TEST_EQ(counter, 0);
            // All are default constructed, distinct instances
            TEST_EQ(prop_type::get(ss).x, -1);
            TEST_EQ(prop_type::get(ss2).x, -1);
            prop_type::get(ss2).x = 2;
            TEST_EQ(prop_type::get(ss).x, -1);
            TEST_EQ(prop_type::get(ss2).x, 2);
            TEST_EQ(counter, 2);
        }
    }
    TEST_EQ(counter, 0);
}
