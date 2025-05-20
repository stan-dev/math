//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

#ifndef BOOST_URL_TEST_RULE_HPP
#define BOOST_URL_TEST_RULE_HPP

#include <boost/url/grammar/charset.hpp>
#include <boost/url/grammar/parse.hpp>
#include <boost/url/grammar/type_traits.hpp>
#include <boost/url/string_view.hpp>
#include <initializer_list>
#include <type_traits>

#include "test_suite.hpp"

namespace boost {
namespace urls {

template<class F>
void
for_each_char(F const& f)
{
    for (int i = 0; i < 256; ++i)
    {
        f(static_cast<char>(i));
    }
}

template<class CharSet>
void
test_char_set(
    CharSet const& cs,
    core::string_view s) noexcept
{
    // each char in s is in the set.
    for(char c : s)
        BOOST_TEST(cs(c));

    // number of chars in
    // set equals s.size()
    std::size_t n = 0;
    for_each_char(
    [&cs, &n](char c)
    {
        if(cs(c))
            ++n;
    });
    BOOST_TEST_EQ(n, s.size());

    // test find_if and find_if_not
    for_each_char(
    [&cs](char c)
    {
        if(cs(static_cast<unsigned char>(0)) || ! cs(c))
        {
            if(cs(c))
            {
                BOOST_TEST(grammar::find_if(
                    &c, &c+1, cs) == &c);
                BOOST_TEST(grammar::find_if_not(
                    &c, &c+1, cs) == &c+1);
            }
            else
            {
                BOOST_TEST(grammar::find_if(
                    &c, &c+1, cs) == &c+1);
                BOOST_TEST(grammar::find_if_not(
                    &c, &c+1, cs) == &c);
            }
            return;
        }

        char buf[40];
        std::memset(
            buf, 0, sizeof(buf));
        buf[19] = c;
        buf[22] = c;
        BOOST_TEST(grammar::find_if(
            buf, buf + sizeof(buf), cs) ==
                buf + 19);
        std::memset(
            buf, c, sizeof(buf));
        buf[19] = 0;
        buf[22] = 0;
        BOOST_TEST(grammar::find_if_not(
            buf, buf + sizeof(buf), cs) ==
                buf + 19);
    });
}

//------------------------------------------------

// rule must match the string
template<class R>
typename std::enable_if<
    grammar::is_rule<R>::value>::type
ok( R const& r,
    core::string_view s)
{
    BOOST_TEST(grammar::parse(s, r).has_value());
}

// rule must match the string and value
template<class R, class V>
typename std::enable_if<
    grammar::is_rule<R>::value>::type
ok( R const& r,
    core::string_view s,
    V const& v)
{
    auto rv = grammar::parse(s, r);
    if(BOOST_TEST(rv.has_value()))
        BOOST_TEST_EQ(rv.value(), v);
}

// rule must fail the string
template<class R>
typename std::enable_if<
    grammar::is_rule<R>::value>::type
bad(
    R const& r,
    core::string_view s)
{
    BOOST_TEST(grammar::parse(s, r).has_error());
}

// rule must fail the string with error
template<class R>
typename std::enable_if<
    grammar::is_rule<R>::value>::type
bad(
    R const& r,
    core::string_view s,
    system::error_code const& e)
{
    auto rv = grammar::parse(s, r);
    if(BOOST_TEST(rv.has_error()))
        BOOST_TEST_EQ(rv.error(), e);
}

} // urls
} // boost

#endif
