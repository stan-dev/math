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

bp::rule<struct seq1_tag, bp::tuple<int, char>> seq1 = "";
bp::rule<struct seq2_tag, bp::tuple<int, char>> seq2 = "";
auto const seq1_def = bp::int_ >> bp::char_('a');
auto const seq2_def = bp::int_ >> bp::char_('b');
BOOST_PARSER_DEFINE_RULES(seq1, seq2);

int main()
{
    // These are just some assorted cases that have, or seemed likely to,
    // cause problems when an internal failure in an alternative wipes out an
    // existing result.  This covers all the cases where "if (!success)"
    // causes the attribute to be overwritten with a default-constructed
    // value.

    {
        auto const parser =
            bp::string("FOO") >> -(bp::string("bar") | bp::string("foo"));

        auto result = bp::parse("FOOfoo", parser);
        BOOST_TEST(result);
        BOOST_TEST(bp::get(*result, bp::llong<0>{}) == std::string("FOO"));
        BOOST_TEST(bp::get(*result, bp::llong<1>{}) == std::string("foo"));
    }

    {
        auto const parser = bp::merge
            [bp::string("FOO") >> (bp::string("bar") | bp::string("foo"))];

        auto result = bp::parse("FOOfoo", parser);
        BOOST_TEST(result);
        BOOST_TEST(*result == std::string("FOOfoo"));
    }

    {
        auto const parser = bp::merge
            [(bp::attr(std::vector<std::string>({"FOO"})) | bp::eps) >>
             (bp::repeat(1)[bp::string("foo")] | bp::eps)];

        auto result = bp::parse("", parser);
        BOOST_TEST(result);
        BOOST_TEST(*result);
        BOOST_TEST(*result == std::vector<std::string>({"FOO"}));
    }

    {
        auto const parser = bp::merge[seq1 >> (seq2 | seq1)];

        auto result = bp::parse("7a9a", parser);
        BOOST_TEST(result);
        BOOST_TEST(*result == (bp::tuple<int, char>(9, 'a')));
    }

    return boost::report_errors();
}
