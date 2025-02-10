//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/json
//

// Test that header file is self-contained.
#include <boost/json/serialize.hpp>

#include <boost/json/parse.hpp>
#include <boost/json/static_resource.hpp>
#include <limits>
#include <sstream>

#include "test_suite.hpp"

namespace boost {
namespace json {

class serialize_test
{
public:
    template<class T>
    static
    std::string
    print(T const& t)
    {
        std::stringstream ss;
        ss << t;
        return ss.str();
    }

    void
    testSerialize()
    {
        {
            value const jv = { 1, 2, 3 };
            BOOST_TEST(serialize(jv) == "[1,2,3]");
            BOOST_TEST(print(jv) == "[1,2,3]");
        }
        {
            array const arr = { 1, 2 ,3 };
            BOOST_TEST(serialize(arr) == "[1,2,3]");
            BOOST_TEST(print(arr) == "[1,2,3]");
        }
        {
            object const obj = { {"k1",1}, {"k2",2} };
            BOOST_TEST(serialize(obj) == "{\"k1\":1,\"k2\":2}");
            BOOST_TEST(print(obj) == "{\"k1\":1,\"k2\":2}");
        }
        {
            string const str = "123";
            BOOST_TEST(serialize(str) == "\"123\"");
            BOOST_TEST(print(str) == "\"123\"");
        }
    }

    void
    testSpecialNumbers()
    {
        using Lims = std::numeric_limits<double>;
        value const jv = {
            Lims::quiet_NaN(), Lims::infinity(), -Lims::infinity() };
        BOOST_TEST( serialize(jv) == "[null,1e99999,-1e99999]" );
        BOOST_TEST( print(jv) == "[null,1e99999,-1e99999]" );

        serialize_options opts;
        opts.allow_infinity_and_nan = true;
        BOOST_TEST( serialize(jv, opts) == "[NaN,Infinity,-Infinity]" );

        std::stringstream ss;
        ss << opts << jv;
        BOOST_TEST( ss.str() == "[NaN,Infinity,-Infinity]" );

        opts.allow_infinity_and_nan = false;
        ss.str("");
        ss << opts << jv;
        BOOST_TEST( ss.str() == "[null,1e99999,-1e99999]" );
    }

    void
    run()
    {
        testSerialize();
        testSpecialNumbers();
    }
};

TEST_SUITE(serialize_test, "boost.json.serialize");

} // namespace json
} // namespace boost
