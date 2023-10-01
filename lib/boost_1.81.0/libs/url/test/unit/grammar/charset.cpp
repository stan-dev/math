//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/grammar/charset.hpp>

#include <boost/url/grammar/alpha_chars.hpp>
#include <boost/url/grammar/parse.hpp>
#include <boost/url/grammar/token_rule.hpp>
#include "test_rule.hpp"

namespace boost {
namespace urls {
namespace grammar {

class charset_test
{
public:
    struct test_chars
    {
        std::size_t& n0;
        std::size_t& n1;

        bool operator()(
            char c) const noexcept
        {
            return c == 'x';
        }

        char const*
        find_if(
            char const* first,
            char const* last) const noexcept
        {
            ++n0;
            while(first != last)
            {
                if(*first == 'x')
                    break;
                ++first;
            }
            return first;
        }

        char const*
        find_if_not(
            char const* first,
            char const* last) const noexcept
        {
            ++n1;
            while(first != last)
            {
                if(*first != 'x')
                    break;
                ++first;
            }
            return first;
        }
    };

    BOOST_STATIC_ASSERT(
        detail::has_find_if<
            test_chars>::value);

    BOOST_STATIC_ASSERT(
        detail::has_find_if_not<
            test_chars>::value);

    void
    testRef()
    {
        BOOST_STATIC_ASSERT(is_charset<
            decltype(ref(alpha_chars))>::value);
        BOOST_TEST(parse("abc", token_rule(
            ref(alpha_chars))).has_value());
    }

    void
    run()
    {
        testRef();

        std::size_t n0 = 0;
        std::size_t n1 = 0;
        test_char_set(
            test_chars{n0, n1}, "x");
        BOOST_TEST_GT(n0, 0u);
        BOOST_TEST_GT(n1, 0u);
    }
};

TEST_SUITE(
    charset_test,
    "boost.url.grammar.charset");

} // grammar
} // urls
} // boost
