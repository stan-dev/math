//
// Copyright (c) 2022 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/grammar/recycled.hpp>

#include "test_suite.hpp"
#include <string>

namespace boost {
namespace urls {
namespace grammar {

struct recycled_test
{
    void
    run()
    {
        // basic
        {
            recycled_ptr<std::string> sp;
            sp->reserve(1000);
            BOOST_TEST(sp->capacity() >= 1000);
        }
        {
            recycled_ptr<std::string> sp;
            BOOST_TEST(sp->capacity() >= 1000);
        }

        // shared
        {
            recycled_ptr<std::string> sp;
            auto sp2 = sp;
            sp->reserve(1000);
            BOOST_TEST(sp->capacity() >= 1000);
            BOOST_TEST(sp2->capacity() >= 1000);
        }
        {
            recycled_ptr<std::string> sp;
            recycled_ptr<std::string> sp2(sp);
            BOOST_TEST(sp->capacity() >= 1000);
            BOOST_TEST(sp2->capacity() >= 1000);
        }

        // coverage
        {
            see_below::recycled_add_impl(1);
            see_below::recycled_remove_impl(1);
        }
    }
};

TEST_SUITE(
    recycled_test,
    "boost.url.grammar.recycled");

} // grammar
} // urls
} // boost
