//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/grammar/string_token.hpp>

#include "test_suite.hpp"

namespace boost {
namespace urls {
namespace grammar {

struct string_token_test
{
    void
    f_impl(
        string_token::arg& dest,
        string_view v)
    {
        char* p = dest.prepare(v.size());
        v.copy(p, v.size());
    }

    template<
        class StringToken = string_token::return_string>
    typename StringToken::result_type
    f(  StringToken&& st = {},
        string_view v = "test")
    {
        f_impl(st, v);
        return st.result();
    }

    void
    run()
    {
        // return_string
        {
            BOOST_TEST_EQ(f(), "test");
        }

        // append_string
        {
            std::string s("url");
            std::string& r(f(string_token::append_to(s)));
            BOOST_TEST_EQ(r, "urltest");
        }

        // assign_string
        {
            std::string s("url");
            std::string& r(f(string_token::assign_to(s)));
            BOOST_TEST_EQ(r, "test");
        }

        // temp_string
        {
            std::string s("url");
            string_view sv;

            sv = f(string_token::preserve_size(s));
            BOOST_TEST_EQ(sv, "test");

            sv = f(string_token::preserve_size(s),
                "supercalifragilisticexpialidocious" );
            BOOST_TEST_EQ(sv,
                "supercalifragilisticexpialidocious");

            sv = f(string_token::preserve_size(s));
            BOOST_TEST_EQ(sv, "test");
        }
    }
};

TEST_SUITE(
    string_token_test,
    "boost.url.grammar.string_token");

} // grammar
} // urls
} // boost
