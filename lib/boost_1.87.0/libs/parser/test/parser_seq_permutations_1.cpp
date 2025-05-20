/**
 *   Copyright (C) 2024 T. Zachary Laine
 *
 *   Distributed under the Boost Software License, Version 1.0. (See
 *   accompanying file LICENSE_1_0.txt or copy at
 *   http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/parser/parser.hpp>

#include <boost/core/lightweight_test.hpp>


namespace bp = boost::parser;

namespace make {
    template<typename... Ts>
    auto tuple(Ts &&... xs)
    {
        return bp::tuple<Ts...>((Ts &&) xs...);
    }
}

/*
{P0, -P0, *P0, P1, P2, P3, eps}
<cartesian product>
{P0, P1, P2, P3, eps, *P0, -P0, (P0 >>/| P2), -(P0 >>/| P1), (-P0 >>/| P1)}

P0 = bp::string("foo");
P1 = bp::string("bar");
P2 = bp::int_;
P3 = bp::char_('c');
*/

using namespace std::literals;

int main()
{
    [[maybe_unused]] int dummy = 0; // for clang-format(!)

    // P0
    {
        auto result =
            bp::parse("foofoo", bp::string("foo") >> bp::string("foo"));
        BOOST_TEST(result);
        BOOST_TEST(*result == (make::tuple("foo"s, "foo"s)));
    }
    {
        auto result =
            bp::parse("foobar", bp::string("foo") >> bp::string("bar"));
        BOOST_TEST(result);
        BOOST_TEST(*result == (make::tuple("foo"s, "bar"s)));
    }
    {
        auto result = bp::parse("foo42", bp::string("foo") >> bp::int_);
        BOOST_TEST(result);
        BOOST_TEST(*result == (make::tuple("foo"s, 42)));
    }
    {
        auto result = bp::parse("fooc", bp::string("foo") >> bp::char_('c'));
        BOOST_TEST(result);
        BOOST_TEST(*result == "fooc"s);
    }
    {
        auto result = bp::parse("foo", bp::string("foo") >> bp::eps);
        BOOST_TEST(result);
        BOOST_TEST(*result == "foo"s);
    }
    {
        auto result =
            bp::parse("foofoofoo", bp::string("foo") >> *bp::string("foo"));
        BOOST_TEST(result);
        BOOST_TEST(*result == std::vector({"foo"s, "foo"s, "foo"s}));
    }
    {
        auto result = bp::parse("foo", bp::string("foo") >> -bp::string("foo"));
        BOOST_TEST(result);
        BOOST_TEST(*result == (make::tuple("foo"s, std::optional<std::string>{})));
    }
    {
        auto result = bp::parse(
            "foofoo42", bp::string("foo") >> (bp::string("foo") >> bp::int_));
        BOOST_TEST(result);
        BOOST_TEST(*result == (make::tuple("foo"s, "foo"s, 42)));
    }
    {
        auto result = bp::parse(
            "foofoobar",
            bp::string("foo") >> -(bp::string("foo") >> bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(
            *result ==
            (make::tuple("foo"s, std::optional(make::tuple("foo"s, "bar"s)))));
    }
    {
        auto result = bp::parse(
            "foofoobar",
            bp::string("foo") >> (-bp::string("foo") >> bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(
            *result == (make::tuple("foo"s, std::optional("foo"s), "bar"s)));
    }
    {
        auto result = bp::parse(
            "foo42", bp::string("foo") >> (bp::string("foo") | bp::int_));
        BOOST_TEST(result);
        BOOST_TEST(
            *result == (make::tuple("foo"s, std::variant<std::string, int>(42))));
    }
    {
        auto result = bp::parse(
            "foo",
            bp::string("foo") >> -(bp::string("foo") | bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(*result == (make::tuple("foo"s, std::optional<std::string>{})));
    }
    {
        auto result = bp::parse(
            "foo",
            bp::string("foo") >> (-bp::string("foo") | bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(
            *result ==
            (make::tuple(
                "foo"s,
                std::variant<std::optional<std::string>, std::string>(
                    std::nullopt))));
    }

    // -P0
    {
        auto result =
            bp::parse("foofoo", -bp::string("foo") >> bp::string("foo"));
        BOOST_TEST(result);
        BOOST_TEST(*result == (make::tuple(std::optional("foo"s), "foo"s)));
    }
    {
        auto result = bp::parse("bar", -bp::string("foo") >> bp::string("bar"));
        BOOST_TEST(result);
        BOOST_TEST(*result == (make::tuple(std::optional<std::string>{}, "bar"s)));
    }
    {
        auto result = bp::parse("42", -bp::string("foo") >> bp::int_);
        BOOST_TEST(result);
        BOOST_TEST(*result == (make::tuple(std::optional<std::string>{}, 42)));
    }
    {
        auto result = bp::parse("c", -bp::string("foo") >> bp::char_('c'));
        BOOST_TEST(result);
        BOOST_TEST(*result == (make::tuple(std::optional<std::string>{}, 'c')));
    }
    {
        auto result = bp::parse("foo", -bp::string("foo") >> bp::eps);
        BOOST_TEST(result);
        BOOST_TEST(*result == std::optional("foo"s));
    }
    {
        auto result =
            bp::parse("foofoo", -bp::string("foo") >> *bp::string("foo"));
        BOOST_TEST(result);
        BOOST_TEST(*result == std::vector({"foo"s, "foo"s}));
    }
    {
        auto result =
            bp::parse("foofoo", -bp::string("foo") >> -bp::string("foo"));
        BOOST_TEST(result);
        BOOST_TEST(
            *result ==
            (make::tuple(std::optional("foo"s), std::optional("foo"s))));
    }
    {
        auto result = bp::parse(
            "foofoo42", -bp::string("foo") >> (bp::string("foo") >> bp::int_));
        BOOST_TEST(result);
        BOOST_TEST(*result == (make::tuple(std::optional("foo"s), "foo"s, 42)));
    }
    {
        auto result = bp::parse(
            "foofoobar",
            -bp::string("foo") >> -(bp::string("foo") >> bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(
            *result ==
            (make::tuple(
                std::optional("foo"s),
                std::optional((make::tuple("foo"s, "bar"s))))));
    }
    {
        auto result = bp::parse(
            "foofoobar",
            -bp::string("foo") >> (-bp::string("foo") >> bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(
            *result ==
            (make::tuple(
                std::optional("foo"s), std::optional("foo"s), "bar"s)));
    }
    {
        auto result = bp::parse(
            "foo42", -bp::string("foo") >> (bp::string("foo") | bp::int_));
        BOOST_TEST(result);
        BOOST_TEST(
            *result ==
            (make::tuple(
                std::optional("foo"s), std::variant<std::string, int>(42))));
    }
    {
        auto result = bp::parse(
            "foobar",
            -bp::string("foo") >> -(bp::string("foo") | bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(
            *result ==
            (make::tuple(std::optional("foo"s), std::optional("bar"s))));
    }
    {
        auto result = bp::parse(
            "foofoo",
            -bp::string("foo") >> (-bp::string("foo") | bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(
            *result ==
            (make::tuple(
                "foo"s,
                std::variant<std::optional<std::string>, std::string>(
                    std::optional("foo"s)))));
    }

    // *P0
    {
        auto result = bp::parse(
            "foo foo",
            bp::lexeme[*bp::string("foo")] >> bp::string("foo"),
            bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(*result == std::vector({"foo"s, "foo"s}));
    }
    {
        auto result =
            bp::parse("foobar", *bp::string("foo") >> bp::string("bar"));
        BOOST_TEST(result);
        BOOST_TEST(*result == std::vector({"foo"s, "bar"s}));
    }
    {
        auto result = bp::parse("foo42", *bp::string("foo") >> bp::int_);
        BOOST_TEST(result);
        BOOST_TEST(*result == (make::tuple(std::vector({"foo"s}), 42)));
    }
    {
        auto result = bp::parse("fooc", *bp::string("foo") >> bp::char_('c'));
        BOOST_TEST(result);
        BOOST_TEST(*result == (make::tuple(std::vector({"foo"s}), 'c')));
    }
    {
        auto result = bp::parse("foo", *bp::string("foo") >> bp::eps);
        BOOST_TEST(result);
        BOOST_TEST(*result == std::vector({"foo"s}));
    }
    {
        auto result = bp::parse(
            "foo foo",
            bp::lexeme[*bp::string("foo")] >> *bp::string("foo"),
            bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(
            *result ==
            (make::tuple(std::vector({"foo"s}), std::vector({"foo"s}))));
    }
    {
        auto result =
            bp::parse("foo", *bp::string("foo") >> -bp::string("foo"));
        BOOST_TEST(result);
        BOOST_TEST(*result == std::vector({"foo"s}));
    }
    {
        auto result = bp::parse(
            "foo foo42",
            bp::lexeme[*bp::string("foo")] >> (bp::string("foo") >> bp::int_),
            bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(*result == (make::tuple(std::vector({"foo"s, "foo"s}), 42)));
    }
    {
        auto result = bp::parse(
            "foo foobar",
            bp::lexeme[*bp::string("foo")] >>
                -(bp::string("foo") >> bp::string("bar")),
            bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(
            *result ==
            (make::tuple(
                std::vector({"foo"s}),
                std::optional((make::tuple("foo"s, "bar"s))))));
    }
    {
        auto result = bp::parse(
            "foo foobar",
            bp::lexeme[*bp::string("foo")] >>
                (-bp::string("foo") >> bp::string("bar")),
            bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(*result == std::vector({"foo"s, "foo"s, "bar"s}));
    }
    {
        auto result = bp::parse(
            "foo 42",
            bp::lexeme[*bp::string("foo")] >> (bp::string("foo") | bp::int_),
            bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(
            *result ==
            (make::tuple(
                std::vector({"foo"s}), std::variant<std::string, int>(42))));
    }
    {
        auto result = bp::parse(
            "foo bar",
            bp::lexeme[*bp::string("foo")] >>
                -(bp::string("foo") | bp::string("bar")),
            bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(*result == std::vector({"foo"s, "bar"s}));
    }
    {
        auto result = bp::parse(
            "foo",
            bp::lexeme[*bp::string("foo")] >>
                (-bp::string("foo") | bp::string("bar")),
            bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(
            *result ==
            (make::tuple(
                std::vector({"foo"s}),
                std::variant<std::optional<std::string>, std::string>(
                    std::nullopt))));
    }

    return boost::report_errors();
}
