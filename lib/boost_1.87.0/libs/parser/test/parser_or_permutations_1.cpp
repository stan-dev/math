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
{P0, P1, P2, P3, eps, *P0, -P0, (P0 |/>> P2), -(P0 |/>> P1), (-P0 |/>> P1)}

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
        auto result = bp::parse("foo", bp::string("foo") | bp::string("foo"));
        BOOST_TEST(result);
        BOOST_TEST(*result == "foo"s);
    }
    {
        auto result = bp::parse("bar", bp::string("foo") | bp::string("bar"));
        BOOST_TEST(result);
        BOOST_TEST(*result == "bar"s);
    }
    {
        auto result = bp::parse("42", bp::string("foo") | bp::int_);
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 1u);
        BOOST_TEST(std::get<int>(*result) == 42);
    }
    {
        auto result = bp::parse("c", bp::string("foo") | bp::char_('c'));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 1u);
        BOOST_TEST(std::get<char>(*result) == 'c');
    }
    {
        auto result = bp::parse("foo", bp::string("foo") | bp::eps);
        BOOST_TEST(result);
        BOOST_TEST(*result == std::optional("foo"s));
    }
    {
        auto result = bp::parse("foo", bp::string("foo") | *bp::string("foo"));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(std::get<std::string>(*result) == "foo"s);
    }
    {
        auto result = bp::parse("", bp::string("foo") | -bp::string("foo"));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 1u);
        BOOST_TEST(std::get<std::optional<std::string>>(*result) == std::nullopt);
    }
    {
        auto result = bp::parse(
            "foo", bp::string("foo") | (bp::string("foo") | bp::int_));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(std::get<std::string>(*result) == "foo"s);
    }
    {
        auto result = bp::parse(
            "bar",
            bp::string("foo") | -(bp::string("foo") | bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 1);
        BOOST_TEST(std::get<std::optional<std::string>>(*result) == "bar"s);
    }
    {
        auto result = bp::parse(
            "", bp::string("foo") | (-bp::string("foo") | bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 1);
        BOOST_TEST(std::get<std::optional<std::string>>(*result) == std::nullopt);
    }
    {
        auto result = bp::parse(
            "foo", bp::string("foo") | (bp::string("foo") >> bp::int_));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(std::get<std::string>(*result) == "foo"s);
    }
    {
        auto result = bp::parse(
            "", bp::string("foo") | -(bp::string("foo") >> bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 1u);
        BOOST_TEST(
            (std::get<std::optional<bp::tuple<std::string, std::string>>>(
                *result)) == std::nullopt);
    }
    {
        auto result = bp::parse(
            "bar",
            bp::string("foo") | (-bp::string("foo") >> bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 1u);
        BOOST_TEST(
            (std::get<bp::tuple<std::optional<std::string>, std::string>>(
                *result)) ==
            (make::tuple(std::optional<std::string>{}, "bar"s)));
    }

    // -P0
    {
        auto result = bp::parse("foo", -bp::string("foo") | bp::string("foo"));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(
            std::get<std::optional<std::string>>(*result) ==
            std::optional("foo"s));
    }
    {
        auto result = bp::parse("", -bp::string("foo") | bp::string("bar"));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(std::get<std::optional<std::string>>(*result) == std::nullopt);
    }
    {
        auto result = bp::parse("", -bp::string("foo") | bp::int_);
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(std::get<std::optional<std::string>>(*result) == std::nullopt);
    }
    {
        auto result = bp::parse("", -bp::string("foo") | bp::char_('c'));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(std::get<std::optional<std::string>>(*result) == std::nullopt);
    }
    {
        auto result = bp::parse("foo", -bp::string("foo") | bp::eps);
        BOOST_TEST(result);
        BOOST_TEST(*result == "foo"s);
    }
    {
        auto result = bp::parse("foo", -bp::string("foo") | *bp::string("foo"));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(std::get<std::optional<std::string>>(*result) == "foo"s);
    }
    {
        auto result = bp::parse("foo", -bp::string("foo") | -bp::string("foo"));
        BOOST_TEST(result);
        BOOST_TEST(*result == "foo"s);
    }
    {
        auto result =
            bp::parse("", -bp::string("foo") | (bp::string("foo") | bp::int_));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(std::get<std::optional<std::string>>(*result) == std::nullopt);
    }
    {
        auto result = bp::parse(
            "foo",
            -bp::string("foo") | -(bp::string("foo") | bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(*result == "foo"s);
    }
    {
        auto result = bp::parse(
            "", -bp::string("foo") | (-bp::string("foo") | bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(std::get<std::optional<std::string>>(*result) == std::nullopt);
    }
    {
        auto result = bp::parse(
            "foo", -bp::string("foo") | (bp::string("foo") >> bp::int_));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(std::get<std::optional<std::string>>(*result) == "foo"s);
    }
    {
        auto result = bp::parse(
            "foo",
            -bp::string("foo") | -(bp::string("foo") >> bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(std::get<std::optional<std::string>>(*result) == "foo"s);
    }
    {
        auto result = bp::parse(
            "foo",
            -bp::string("foo") | (-bp::string("foo") >> bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(std::get<std::optional<std::string>>(*result) == "foo"s);
    }

    // *P0
    {
        auto result = bp::parse(
            "foo", bp::lexeme[*bp::string("foo")] | bp::string("foo"));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(
            std::get<std::vector<std::string>>(*result) == std::vector({"foo"s}));
    }
    {
        auto result =
            bp::parse("foofoo", *bp::string("foo") | bp::string("bar"));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(
            std::get<std::vector<std::string>>(*result) ==
            std::vector({"foo"s, "foo"s}));
    }
    {
        auto result = bp::parse("", *bp::string("foo") | bp::int_);
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(
            std::get<std::vector<std::string>>(*result) ==
            std::vector<std::string>{});
    }
    {
        auto result = bp::parse("", *bp::string("foo") | bp::char_('c'));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(
            std::get<std::vector<std::string>>(*result) ==
            std::vector<std::string>{});
    }
    {
        auto result = bp::parse("foofoo", *bp::string("foo") | bp::eps);
        BOOST_TEST(result);
        BOOST_TEST(*result == std::optional(std::vector({"foo"s, "foo"s})));
    }
    {
        auto result = bp::parse(
            "foofoo", bp::lexeme[*bp::string("foo")] | *bp::string("foo"));
        BOOST_TEST(result);
        BOOST_TEST(*result == std::vector({"foo"s, "foo"s}));
    }
    {
        auto result =
            bp::parse("foofoo", *bp::string("foo") | -bp::string("foo"));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(
            std::get<std::vector<std::string>>(*result) ==
            std::vector({"foo"s, "foo"s}));
    }
    {
        auto result = bp::parse(
            "foofoo",
            bp::lexeme[*bp::string("foo")] | (bp::string("foo") | bp::int_));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(
            std::get<std::vector<std::string>>(*result) ==
            std::vector({"foo"s, "foo"s}));
    }
    {
        auto result = bp::parse(
            "foofoo",
            bp::lexeme[*bp::string("foo")] |
                -(bp::string("foo") | bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(
            std::get<std::vector<std::string>>(*result) ==
            std::vector({"foo"s, "foo"s}));
    }
    {
        auto result = bp::parse(
            "foofoo",
            bp::lexeme[*bp::string("foo")] |
                (-bp::string("foo") | bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(
            std::get<std::vector<std::string>>(*result) ==
            std::vector({"foo"s, "foo"s}));
    }
    {
        auto result = bp::parse(
            "foofoo",
            bp::lexeme[*bp::string("foo")] | (bp::string("foo") >> bp::int_));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(
            std::get<std::vector<std::string>>(*result) ==
            std::vector({"foo"s, "foo"s}));
    }
    {
        auto result = bp::parse(
            "foofoo",
            bp::lexeme[*bp::string("foo")] |
                -(bp::string("foo") >> bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(
            std::get<std::vector<std::string>>(*result) ==
            std::vector({"foo"s, "foo"s}));
    }
    {
        auto result = bp::parse(
            "foofoo",
            bp::lexeme[*bp::string("foo")] |
                (-bp::string("foo") >> bp::string("bar")));
        BOOST_TEST(result);
        BOOST_TEST(result->index() == 0);
        BOOST_TEST(
            std::get<std::vector<std::string>>(*result) ==
            std::vector({"foo"s, "foo"s}));
    }

return boost::report_errors();
}
