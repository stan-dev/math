//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/grammar/error.hpp>

#include "test_suite.hpp"

#include <memory>

namespace boost {
namespace urls {
namespace grammar {

class error_test
{
public:
    void check(error e)
    {
        auto const ec = BOOST_URL_ERR(e);
        BOOST_TEST_NE(ec.category().name(), nullptr);
        BOOST_TEST(! ec.message().empty());
        BOOST_TEST(ec.category().default_error_condition(
            static_cast<int>(e)).category() == ec.category());
    }

    void check(condition c, error e)
    {
        {
            auto const ec = BOOST_URL_ERR(e);
            BOOST_TEST_NE(ec.category().name(), nullptr);
            BOOST_TEST(! ec.message().empty());
            BOOST_TEST_EQ(ec, c);
        }
        {
            auto ec = make_error_condition(c);
            BOOST_TEST_NE(ec.category().name(), nullptr);
            BOOST_TEST(! ec.message().empty());
            BOOST_TEST_EQ(ec, c);
        }
    }

    void
    run()
    {
        check(error::need_more);
        check(error::mismatch);
        check(error::end_of_range);
        check(error::leftover);

        check(condition::fatal, error::invalid);
        check(condition::fatal, error::out_of_range);
    }
};

TEST_SUITE(
    error_test,
    "boost.url.grammar.error");

} // grammar
} // urls
} // boost
