//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/json
//

// Test that header file is self-contained.
#include <boost/json/value_to.hpp>

#include "test_suite.hpp"

#include <map>
#include <unordered_map>
#include <vector>

BOOST_JSON_NS_BEGIN

class value_to_test
{
public:
    
    template<class T>
    void
    check(T t)
    {
        BOOST_TEST(value_to<T>(value_from(t)) == t);
    }

    void
    testNumberCast()
    {
        check((short)-1);
        check((int)-2);
        check((long)-3);
        check((long long)-4);
        check((unsigned short)1);
        check((unsigned int)2);
        check((unsigned long)3);
        check((unsigned long long)4);
        check((float)1.5);
        check((double)2.5);
        check(true);
    }

    void
    testJsonTypes()
    {
        value_to<object>(value(object_kind));
        value_to<array>(value(array_kind));
        value_to<string>(value(string_kind));
    }

    void
    testGenerics()
    {
        check(std::string("test"));
        check(std::map<std::string, int>
        {
            {"a", 1}, {"b", 2}, {"c", 3}
        });
        check(std::unordered_map<std::string, int>
        {
            { "a", 1 }, {"b", 2}, {"c", 3}
        });
        check(std::vector<int>{1, 2, 3, 4});
    }

    void
    run()
    {
        testNumberCast();
        testJsonTypes();
        testGenerics();
    }
};

TEST_SUITE(value_to_test, "boost.json.value_to");

BOOST_JSON_NS_END
