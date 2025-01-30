/**
 *   Copyright (C) 2018 T. Zachary Laine
 *
 *   Distributed under the Boost Software License, Version 1.0. (See
 *   accompanying file LICENSE_1_0.txt or copy at
 *   http://www.boost.org/LICENSE_1_0.txt)
 */
#include <boost/parser/parser.hpp>
#include <boost/parser/transcode_view.hpp>

#include "ill_formed.hpp"

#include <boost/core/lightweight_test.hpp>

#include <any>
#include <deque>


using namespace boost::parser;


constexpr callback_rule<struct callback_char_rule_tag, char>
    callback_char_rule = "callback_char_rule";
constexpr auto callback_char_rule_def = char_;
BOOST_PARSER_DEFINE_RULES(callback_char_rule);

struct callback_char_rule_tag
{};

int main()
{

// full_parse_api
{
    std::string const str = "a";

    // attr out param, iter/sent
    {
        char out = 0;
        auto first = str.c_str();
        BOOST_TEST(prefix_parse(
            first, boost::parser::detail::text::null_sentinel, char_, out));
        first = str.c_str();
        BOOST_TEST(out == 'a');
        out = 0;
        first = str.c_str();
        BOOST_TEST(!prefix_parse(
            first,
            boost::parser::detail::text::null_sentinel,
            char_('b'),
            out));
        BOOST_TEST(out == 0);
    }
    // attr out param, range
    {
        char out = 0;
        BOOST_TEST(parse(str, char_, out));
        BOOST_TEST(out == 'a');
        out = 0;
        BOOST_TEST(!parse(str, char_('b'), out));
        BOOST_TEST(out == 0);
    }
    // attr out param, pointer-as-range
    {
        char out = 0;
        BOOST_TEST(parse(null_term(str.c_str()), char_, out));
        BOOST_TEST(out == 'a');
        out = 0;
        BOOST_TEST(!parse(null_term(str.c_str()), char_('b'), out));
        BOOST_TEST(out == 0);
    }

    // returned attr, iter/sent
    {
        auto first = str.c_str();
        BOOST_TEST(prefix_parse(
            first, boost::parser::detail::text::null_sentinel, char_));
        first = str.c_str();
        BOOST_TEST(
            *prefix_parse(
                first, boost::parser::detail::text::null_sentinel, char_) ==
            'a');
        first = str.c_str();
        BOOST_TEST(!prefix_parse(
            first, boost::parser::detail::text::null_sentinel, char_('b')));
    }
    // returned attr, range
    {
        BOOST_TEST(parse(str, char_));
        BOOST_TEST(*parse(str, char_) == 'a');
        BOOST_TEST(!parse(str, char_('b')));
    }
    // returned attr, pointer-as-range
    {
        BOOST_TEST(parse(null_term(str.c_str()), char_));
        BOOST_TEST(*parse(null_term(str.c_str()), char_) == 'a');
        BOOST_TEST(!parse(null_term(str.c_str()), char_('b')));
    }
    // returned attr, UTF-16
    {
        BOOST_TEST(parse(u"a", char_));
        auto const result = *parse(u"a", char_);
        BOOST_TEST(uint16_t(result) == uint16_t('a'));
        BOOST_TEST(!parse(u"a", char_('b')));
    }

    // attr out param, using skipper, iter/sent
    {
        char out = 0;
        auto first = str.c_str();
        BOOST_TEST(prefix_parse(
            first,
            boost::parser::detail::text::null_sentinel,
            char_,
            ws,
            out));
        first = str.c_str();
        BOOST_TEST(out == 'a');
        out = 0;
        first = str.c_str();
        auto ws_copy = ws;
        BOOST_TEST(!prefix_parse(
            first,
            boost::parser::detail::text::null_sentinel,
            char_('b'),
            ws_copy,
            out));
        BOOST_TEST(out == 0);
    }
    // attr out param, using skipper, range
    {
        char out = 0;
        BOOST_TEST(parse(str, char_, ws, out));
        BOOST_TEST(out == 'a');
        out = 0;
        BOOST_TEST(!parse(str, char_('b'), ws, out));
        BOOST_TEST(out == 0);
    }
    // attr out param, using skipper, pointer-as-range
    {
        char out = 0;
        BOOST_TEST(parse(null_term(str.c_str()), char_, ws, out));
        BOOST_TEST(out == 'a');
        out = 0;
        auto ws_copy = ws;
        BOOST_TEST(!parse(null_term(str.c_str()), char_('b'), ws_copy, out));
        BOOST_TEST(out == 0);
    }

    // returned attr, using skipper, iter/sent
    {
        auto first = str.c_str();
        BOOST_TEST(prefix_parse(
            first, boost::parser::detail::text::null_sentinel, char_, ws));
        first = str.c_str();
        BOOST_TEST(
            *prefix_parse(
                first, boost::parser::detail::text::null_sentinel, char_, ws) ==
            'a');
        first = str.c_str();
        auto ws_copy = ws;
        BOOST_TEST(!prefix_parse(
            first,
            boost::parser::detail::text::null_sentinel,
            char_('b'),
            ws_copy));
    }
    // returned attr, using skipper, range
    {
        BOOST_TEST(parse(str, char_, ws));
        BOOST_TEST(*parse(str, char_, ws) == 'a');
        BOOST_TEST(!parse(str, char_('b'), ws));
    }
    // returned attr, using skipper, pointer-as-range
    {
        BOOST_TEST(parse(null_term(str.c_str()), char_, ws));
        BOOST_TEST(*parse(null_term(str.c_str()), char_, ws) == 'a');
        auto ws_copy = ws;
        BOOST_TEST(!parse(null_term(str.c_str()), char_('b'), ws_copy));
    }

    // callback, iter/sent
    {
        char out = 0;
        auto callbacks = [&out](auto tag, auto x) { out = x; };
        auto first = str.c_str();
        BOOST_TEST(callback_prefix_parse(
            first,
            boost::parser::detail::text::null_sentinel,
            callback_char_rule,
            callbacks));
        first = str.c_str();
        BOOST_TEST(out == 'a');
    }
    // callback, range
    {
        char out = 0;
        auto callbacks = [&out](auto tag, auto x) { out = x; };
        BOOST_TEST(callback_parse(str, callback_char_rule, callbacks));
        BOOST_TEST(out == 'a');
    }
    // callback, pointer-as-range
    {
        char out = 0;
        auto callbacks = [&out](auto tag, auto x) { out = x; };
        BOOST_TEST(callback_parse(null_term(str.c_str()), callback_char_rule, callbacks));
        BOOST_TEST(out == 'a');
    }

    // callback, using skipper, iter/sent
    {
        char out = 0;
        auto callbacks = [&out](auto tag, auto x) { out = x; };
        auto first = str.c_str();
        BOOST_TEST(callback_prefix_parse(
            first,
            boost::parser::detail::text::null_sentinel,
            callback_char_rule,
            ws,
            callbacks));
        first = str.c_str();
        BOOST_TEST(out == 'a');
    }
    {
        char out = 0;
        auto callbacks = [&out](auto tag, auto x) { out = x; };
        auto first = str.c_str();
        auto ws_copy = ws;
        BOOST_TEST(callback_prefix_parse(
            first,
            boost::parser::detail::text::null_sentinel,
            callback_char_rule,
            ws_copy,
            callbacks));
        first = str.c_str();
        BOOST_TEST(out == 'a');
    }
    // callback, using skipper, range
    {
        char out = 0;
        auto callbacks = [&out](auto tag, auto x) { out = x; };
        BOOST_TEST(
            callback_parse(str, callback_char_rule, ws, callbacks));
        BOOST_TEST(out == 'a');
    }
    // callback, using skipper, pointer-as-range
    {
        char out = 0;
        auto callbacks = [&out](auto tag, auto x) { out = x; };
        auto ws_copy = ws;
        BOOST_TEST(callback_parse(
            null_term(str.c_str()), callback_char_rule, ws_copy, callbacks));
        BOOST_TEST(out == 'a');
    }
}

// basic
{
    constexpr auto parser_1 = char_ >> char_;
    constexpr auto parser_2 = char_ >> char_ >> char_;
    constexpr auto parser_3 = char_ | char_;
    constexpr auto parser_4 = char_('a') | char_('b') | char_('c');
    constexpr auto parser_5 = char_('a') | char_('b') | eps;

    {
        char const * str = "a";
        BOOST_TEST(parse(null_term(str), char_));
        BOOST_TEST(!parse(null_term(str), char_('b')));
    }
    {
        char const * str = "a";
        char c = '\0';
        BOOST_TEST(parse(null_term(str), char_, c));
        BOOST_TEST(c == 'a');
        BOOST_TEST(!parse(null_term(str), char_('b')));
    }
    {
        char const * str = "b";
        char c = '\0';
        BOOST_TEST(parse(null_term(str), char_("ab"), c));
        BOOST_TEST(c == 'b');
        BOOST_TEST(!parse(null_term(str), char_("cd")));
    }
    {
        char const * str = "b";
        char c = '\0';
        std::string const pattern_1 = "ab";
        std::string const pattern_2 = "cd";
        BOOST_TEST(parse(null_term(str), char_(pattern_1), c));
        BOOST_TEST(c == 'b');
        BOOST_TEST(!parse(null_term(str), char_(pattern_2)));
    }
    {
        char const * str = "b";
        char c = '\0';
        BOOST_TEST(parse(null_term(str), char_('a', 'b'), c));
        BOOST_TEST(c == 'b');
        BOOST_TEST(!parse(null_term(str), char_('c', 'd')));
    }
    {
        char const * str = " ";
        BOOST_TEST(parse(null_term(str), blank));
        BOOST_TEST(!parse(null_term(str), lower));
    }
    {
        char const * str = "ab";
        BOOST_TEST(!parse(null_term(str), char_));
        BOOST_TEST(parse(null_term(str), parser_1));
        BOOST_TEST(!parse(null_term(str), parser_2));
    }
    {
        std::string str = "ab";
        auto first = str.c_str();
        BOOST_TEST(prefix_parse(
            first, boost::parser::detail::text::null_sentinel, char_));
        first = str.c_str();
        BOOST_TEST(prefix_parse(
            first, boost::parser::detail::text::null_sentinel, parser_1));
        first = str.c_str();
        BOOST_TEST(!parse(str, parser_2));
    }
    {
        char const * str = "ab";
        tuple<char, char> result;
        BOOST_TEST(parse(null_term(str), parser_1, result));
        using namespace boost::parser::literals;
        BOOST_TEST(get(result, 0_c) == 'b');
        BOOST_TEST(get(result, 1_c) == '\0');
    }
    {
        char const * str = "abc";
        BOOST_TEST(!parse(null_term(str), parser_1));
        BOOST_TEST(parse(null_term(str), parser_2));
    }
    {
        std::string str = "abc";
        auto first = str.c_str();
        BOOST_TEST(prefix_parse(
            first, boost::parser::detail::text::null_sentinel, parser_1));
        first = str.c_str();
        BOOST_TEST(prefix_parse(
            first, boost::parser::detail::text::null_sentinel, parser_2));
    }
    {
        char const * str = "abc";
        tuple<char, char, char> result;
        BOOST_TEST(parse(null_term(str), parser_2, result));
        using namespace boost::parser::literals;
        BOOST_TEST(get(result, 0_c) == 'c');
        BOOST_TEST(get(result, 1_c) == '\0');
        BOOST_TEST(get(result, 2_c) == '\0');
    }
    {
        char const * str = "a";
        BOOST_TEST(parse(null_term(str), parser_3));
        BOOST_TEST(parse(null_term(str), parser_4));
    }
    {
        char const * str = "a";
        char c = '\0';
        BOOST_TEST(parse(null_term(str), parser_3, c));
        BOOST_TEST(c == 'a');
    }
    {
        char const * str = "a";
        char c = '\0';
        BOOST_TEST(parse(null_term(str), parser_4, c));
        BOOST_TEST(c == 'a');
    }
    {
        char const * str = "z";
        BOOST_TEST(parse(null_term(str), parser_3));
        BOOST_TEST(!parse(null_term(str), parser_4));
    }
    {
        char const * str = "a";
        BOOST_TEST(parse(null_term(str), parser_5));
    }
    {
        char const * str = "z";
        BOOST_TEST(!parse(null_term(str), parser_5));
    }
    {
        std::string str = "z";
        auto first = str.c_str();
        BOOST_TEST(prefix_parse(
            first, boost::parser::detail::text::null_sentinel, parser_5));
    }
    {
        char const * str = "a";
        std::optional<char> c;
        BOOST_TEST(parse(null_term(str), parser_5, c));
        BOOST_TEST(c == 'a');
    }
    {
        char const * str = "z";
        std::optional<char> c;
        BOOST_TEST(!parse(null_term(str), parser_5, c));
    }
    {
        std::string str = "z";
        std::optional<char> c;
        auto first = str.c_str();
        BOOST_TEST(prefix_parse(
            first, boost::parser::detail::text::null_sentinel, parser_5, c));
        BOOST_TEST(c == std::nullopt);
    }
}

// int_
{
    {
        char const * str = "-42";
        short i = 0;
        BOOST_TEST(parse(null_term(str), short_, i));
        BOOST_TEST(i == -42);
    }
    {
        char const * str = "42";
        short i = 0;
        BOOST_TEST(parse(null_term(str), short_, i));
        BOOST_TEST(i == 42);
    }
    {
        char const * str = "-2000000000";
        int i = 0;
        BOOST_TEST(parse(null_term(str), int_, i));
        BOOST_TEST(i == -2000000000);
    }
    {
        char const * str = "2000000000";
        int i = 0;
        BOOST_TEST(parse(null_term(str), int_, i));
        BOOST_TEST(i == 2000000000);
    }
    {
        char const * str = "-2000000000";
        long i = 0;
        BOOST_TEST(parse(null_term(str), long_, i));
        BOOST_TEST(i == -2000000000);
    }
    {
        char const * str = "2000000000";
        long i = 0;
        BOOST_TEST(parse(null_term(str), long_, i));
        BOOST_TEST(i == 2000000000);
    }
    {
        char const * str = "-4000000000";
        long long i = 0;
        BOOST_TEST(parse(null_term(str), long_long, i));
        BOOST_TEST(i == -4000000000LL);
    }
    {
        char const * str = "4000000000";
        long long i = 0;
        BOOST_TEST(parse(null_term(str), long_long, i));
        BOOST_TEST(i == 4000000000LL);
    }
}

// uint_
{
    {
        char const * str = "10011";
        unsigned int i = 0;
        BOOST_TEST(parse(null_term(str), bin, i));
        BOOST_TEST(i == 19);
    }
    {
        char const * str = "107";
        unsigned int i = 0;
        BOOST_TEST(parse(null_term(str), oct, i));
        BOOST_TEST(i == 71);
    }
    {
        char const * str = "beef";
        unsigned int i = 0;
        BOOST_TEST(parse(null_term(str), hex, i));
        BOOST_TEST(i == 48879);
    }

    {
        char const * str = "42";
        unsigned int i = 0;
        BOOST_TEST(parse(null_term(str), ushort_, i));
        BOOST_TEST(i == 42);
    }
    {
        char const * str = "-42";
        unsigned int i = 3;
        BOOST_TEST(!parse(null_term(str), uint_, i));
        BOOST_TEST(i == 0);
    }
    {
        char const * str = "42";
        unsigned int i = 0;
        BOOST_TEST(parse(null_term(str), uint_, i));
        BOOST_TEST(i == 42);
    }
    {
        char const * str = "42";
        unsigned long i = 0;
        BOOST_TEST(parse(null_term(str), ulong_, i));
        BOOST_TEST(i == 42);
    }
    {
        char const * str = "42";
        unsigned long long i = 0;
        BOOST_TEST(parse(null_term(str), ulong_long, i));
        BOOST_TEST(i == 42);
    }
}

// bool_
{
    {
        char const * str = "";
        bool b = false;
        BOOST_TEST(!parse(null_term(str), bool_, b));
    }
    {
        char const * str = "true";
        bool b = false;
        BOOST_TEST(parse(null_term(str), bool_, b));
        BOOST_TEST(b == true);
    }
    {
        char const * str = "false ";
        bool b = true;
        BOOST_TEST(!parse(null_term(str), bool_, b));
    }
    {
        std::string str = "false ";
        bool b = true;
        auto first = str.c_str();
        BOOST_TEST(prefix_parse(
            first, boost::parser::detail::text::null_sentinel, bool_, b));
        BOOST_TEST(b == false);
    }
    {
        char const * str = "true ";
        auto r = boost::parser::as_utf32(boost::parser::null_term(str));
        bool b = false;
        auto first = r.begin();
        auto const last = r.end();
        BOOST_TEST(prefix_parse(first, last, bool_, b));
        BOOST_TEST(b == true);
    }
    {
        char const * str = "false";
        auto r = boost::parser::as_utf32(boost::parser::null_term(str));
        bool b = true;
        auto first = r.begin();
        auto const last = r.end();
        BOOST_TEST(prefix_parse(first, last, bool_, b));
        BOOST_TEST(b == false);
    }
}

// star
{
    {
        constexpr auto parser = *char_;
        {
            char const * str = "";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>());
        }
        {
            char const * str = "a";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>({'a'}));
        }
        {
            char const * str = "ba";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>({'b', 'a'}));
        }
    }

    {
        constexpr auto parser = *char_('b');
        {
            char const * str = "";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>());
        }
        {
            char const * str = "b";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>({'b'}));
        }
        {
            char const * str = "bb";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>({'b', 'b'}));
        }
    }
}

// plus
{
    {
        constexpr auto parser = +char_;

        {
            char const * str = "";
            std::vector<char> chars;
            BOOST_TEST(!parse(null_term(str), parser, chars));
        }
        {
            char const * str = "a";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>({'a'}));
        }
        {
            char const * str = "ba";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>({'b', 'a'}));
        }
    }

    {
        constexpr auto parser = +char_('b');

        {
            char const * str = "";
            std::vector<char> chars;
            BOOST_TEST(!parse(null_term(str), parser, chars));
        }
        {
            char const * str = "b";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>({'b'}));
        }
        {
            char const * str = "bb";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>({'b', 'b'}));
        }
    }
}

// star_and_plus_collapsing
{
    {
        constexpr auto parser = +(+char_('b'));

        {
            char const * str = "";
            std::vector<char> chars;
            BOOST_TEST(!parse(null_term(str), parser, chars));
        }
        {
            char const * str = "b";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>({'b'}));
        }
        {
            char const * str = "bb";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>({'b', 'b'}));
        }
    }

    {
        constexpr auto parser = **char_('z');

        {
            char const * str = "";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>());
        }
        {
            char const * str = "z";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>({'z'}));
        }
        {
            char const * str = "zz";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>({'z', 'z'}));
        }
    }

    {
        constexpr auto parser = +*char_('z');

        {
            char const * str = "";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>());
        }
        {
            char const * str = "z";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>({'z'}));
        }
        {
            char const * str = "zz";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>({'z', 'z'}));
        }
    }

    {
        constexpr auto parser = *+char_('z');

        {
            char const * str = "";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>());
        }
        {
            char const * str = "z";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>({'z'}));
        }
        {
            char const * str = "zz";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>({'z', 'z'}));
        }
    }
}

// action
{
    {{char const * str = "";
    std::stringstream ss;
    auto action = [&ss](auto & context) { ss << _attr(context); };
    auto parser = *char_('b')[action];
    BOOST_TEST(parse(null_term(str), parser));
    BOOST_TEST(ss.str() == "");
}
{
    char const * str = "b";
    std::stringstream ss;
    auto action = [&ss](auto & context) { ss << _attr(context); };
    auto parser = *char_('b')[action];
    BOOST_TEST(parse(null_term(str), parser));
    BOOST_TEST(ss.str() == "b");
}
{
    char const * str = "bb";
    std::stringstream ss;
    auto action = [&ss](auto & context) { ss << _attr(context); };
    auto parser = *char_('b')[action];
    BOOST_TEST(parse(null_term(str), parser));
    BOOST_TEST(parse(null_term(str), parser));
    BOOST_TEST(ss.str() == "bbbb");
}
}

{
    {
        char const * str = "";
        std::stringstream ss;
        auto action = [&ss](auto & context) { ss << _attr(context); };
        auto parser = +char_('b')[action];
        BOOST_TEST(!parse(null_term(str), parser));
        BOOST_TEST(ss.str() == "");
    }
    {
        char const * str = "b";
        std::stringstream ss;
        auto action = [&ss](auto & context) { ss << _attr(context); };
        auto parser = +char_('b')[action];
        BOOST_TEST(parse(null_term(str), parser));
        BOOST_TEST(ss.str() == "b");
    }
    {
        char const * str = "bb";
        std::stringstream ss;
        auto action = [&ss](auto & context) { ss << _attr(context); };
        auto parser = +char_('b')[action];
        BOOST_TEST(parse(null_term(str), parser));
        BOOST_TEST(parse(null_term(str), parser));
        BOOST_TEST(ss.str() == "bbbb");
    }
}
}

// star_as_string_or_vector
{
    {
        constexpr auto parser = *char_('z');

        {
            char const * str = "";
            std::string chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == "");
        }
        {
            char const * str = "z";
            std::string chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == "z");
        }
        {
            char const * str = "zz";
            std::string chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == "zz");
        }
    }

    {
        constexpr auto parser = *char_('z');

        {
            char const * str = "";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>());
        }
        {
            char const * str = "z";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>({'z'}));
        }
        {
            char const * str = "zz";
            std::vector<char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<char>({'z', 'z'}));
        }
    }

    {
        constexpr auto parser = *string("zs");

        {
            char const * str = "";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<std::string>{});

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(null_term(str), parser);
                BOOST_TEST(chars);
                BOOST_TEST(chars->empty());
            }
        }
        {
            std::string str = "z";
            {
                std::vector<std::string> chars;
                BOOST_TEST(!parse(str, parser, chars));
            }
            {
                std::vector<std::string> chars;
                auto first = str.c_str();
                BOOST_TEST(prefix_parse(
                    first,
                    boost::parser::detail::text::null_sentinel,
                    parser,
                    chars));
                BOOST_TEST(chars == std::vector<std::string>{});
            }

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser);
                BOOST_TEST(!chars);
            }
            {
                auto first = str.c_str();
                std::optional<std::vector<std::string>> const chars =
                    prefix_parse(
                        first,
                        boost::parser::detail::text::null_sentinel,
                        parser);
                BOOST_TEST(chars);
                BOOST_TEST(chars->empty());
            }
        }
        {
            char const * str = "zs";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<std::string>({"zs"}));

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(null_term(str), parser);
                BOOST_TEST(chars);
                BOOST_TEST(*chars == std::vector<std::string>({"zs"}));
            }
        }
        {
            char const * str = "zszs";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<std::string>({"zs", "zs"}));

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(null_term(str), parser);
                BOOST_TEST(chars);
                BOOST_TEST(*chars == std::vector<std::string>({"zs", "zs"}));
            }
        }
    }

    {
        constexpr auto parser = *string("zs");

        {
            char const * str = "";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<std::string>{});
        }
        {
            char const * str = "z";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(null_term(str), parser, chars));
        }
        {
            std::string str = "z";
            std::vector<std::string> chars;
            auto first = str.c_str();
            BOOST_TEST(prefix_parse(
                first,
                boost::parser::detail::text::null_sentinel,
                parser,
                chars));
        }
        {
            char const * str = "zs";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<std::string>({"zs"}));
        }
        {
            char const * str = "zszs";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<std::string>({"zs", "zs"}));
        }
    }
}

// omit
{
    {
        constexpr auto parser = omit[*+char_('z')];

        {
            char const * str = "";
            BOOST_TEST(parse(null_term(str), parser));
        }
        {
            char const * str = "z";
            BOOST_TEST(parse(null_term(str), parser));
        }
        {
            char const * str = "zz";
            BOOST_TEST(parse(null_term(str), parser));
        }
        {
            char const * str = "";
            BOOST_TEST(parse(null_term(str), parser));
        }
        {
            char const * str = "z";
            BOOST_TEST(parse(null_term(str), parser));
        }
        {
            char const * str = "zz";
            BOOST_TEST(parse(null_term(str), parser));
        }
    }

    {
        constexpr auto parser = omit[*string("zs")];

        {
            char const * str = "";
            BOOST_TEST(parse(null_term(str), parser));
        }
        {
            char const * str = "z";
            BOOST_TEST(!parse(null_term(str), parser));
        }
        {
            std::string str = "z";
            auto first = str.c_str();
            BOOST_TEST(prefix_parse(
                first, boost::parser::detail::text::null_sentinel, parser));
        }
        {
            char const * str = "zs";
            BOOST_TEST(parse(null_term(str), parser));
        }
        {
            char const * str = "zszs";
            BOOST_TEST(parse(null_term(str), parser));
        }
    }
}

// repeat
{
    {
        constexpr auto parser = repeat(2, 3)[string("zs")];

        {
            char const * str = "";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<std::string>{});

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(null_term(str), parser);
                BOOST_TEST(!chars);
            }
        }
        {
            char const * str = "z";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<std::string>{});

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(null_term(str), parser);
                BOOST_TEST(!chars);
            }
        }
        {
            char const * str = "zs";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<std::string>{});

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(null_term(str), parser);
                BOOST_TEST(!chars);
            }
        }
        {
            char const * str = "zszs";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<std::string>({"zs", "zs"}));

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(null_term(str), parser);
                BOOST_TEST(chars);
                BOOST_TEST(*chars == std::vector<std::string>({"zs", "zs"}));
            }
        }
    }
}

// raw
{
    {
        constexpr auto parser = raw[*string("zs")];
        using range_t = BOOST_PARSER_SUBRANGE<std::string::const_iterator>;

        {
            std::string const str = "";
            range_t r;
            BOOST_TEST(parse(str, parser, r));
            BOOST_TEST(r.begin() == str.begin());
            BOOST_TEST(r.end() == str.begin());
        }
        {
            std::string const str = "z";
            range_t r;
            BOOST_TEST(!parse(str, parser, r));
        }
        {
            std::string const str = "z";
            range_t r;
            auto first = str.begin();
            BOOST_TEST(prefix_parse(first, str.end(), parser, r));
            BOOST_TEST(r.begin() == str.begin());
            BOOST_TEST(r.end() == str.begin());
        }
        {
            std::string const str = "zs";
            range_t r;
            BOOST_TEST(parse(str, parser, r));
            BOOST_TEST(r.begin() == str.begin());
            BOOST_TEST(r.end() == str.end());
        }
        {
            std::string const str = "zszs";
            range_t r;
            BOOST_TEST(parse(str, parser, r));
            BOOST_TEST(r.begin() == str.begin());
            BOOST_TEST(r.end() == str.end());
        }
        {
            std::string const str = "";
            std::optional<range_t> result = parse(str, parser);
            BOOST_TEST(result);
            BOOST_TEST(result->begin() == str.begin());
            BOOST_TEST(result->end() == str.begin());
        }
        {
            std::string const str = "z";
            std::optional<range_t> result = parse(str, parser);
            BOOST_TEST(!result);
        }
        {
            std::string const str = "z";
            auto first = str.begin();
            std::optional<range_t> result =
                prefix_parse(first, str.end(), parser);
            BOOST_TEST(result);
            BOOST_TEST(result->begin() == str.begin());
            BOOST_TEST(result->end() == str.begin());
        }
        {
            std::string const str = "zs";
            std::optional<range_t> result = parse(str, parser);
            BOOST_TEST(result);
            BOOST_TEST(result->begin() == str.begin());
            BOOST_TEST(result->end() == str.end());
        }
        {
            std::string const str = "zszs";
            std::optional<range_t> result = parse(str, parser);
            BOOST_TEST(result);
            BOOST_TEST(result->begin() == str.begin());
            BOOST_TEST(result->end() == str.end());
        }
    }
}

// delimited
{
    {
        constexpr auto parser = string("yay") % ',';

        {
            char const * str = "";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<std::string>{});

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(null_term(str), parser);
                BOOST_TEST(!chars);
            }
        }
        {
            char const * str = "z";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<std::string>{});

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(null_term(str), parser);
                BOOST_TEST(!chars);
            }
        }
        {
            char const * str = ",";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<std::string>{});

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(null_term(str), parser);
                BOOST_TEST(!chars);
            }
        }
        {
            char const * str = ",yay";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<std::string>{});

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(null_term(str), parser);
                BOOST_TEST(!chars);
            }
        }
        {
            char const * str = "yay";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<std::string>({"yay"}));

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(null_term(str), parser);
                BOOST_TEST(chars);
                BOOST_TEST(*chars == std::vector<std::string>({"yay"}));
            }
        }
        {
            std::string str = "yayyay";
            {
                std::vector<std::string> chars;
                BOOST_TEST(!parse(str, parser, chars));
            }
            {
                std::vector<std::string> chars;
                auto first = str.c_str();
                BOOST_TEST(prefix_parse(
                    first,
                    boost::parser::detail::text::null_sentinel,
                    parser,
                    chars));
                BOOST_TEST(chars == std::vector<std::string>({"yay"}));
            }

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser);
                BOOST_TEST(!chars);
            }
            {
                auto first = str.c_str();
                std::optional<std::vector<std::string>> const chars =
                    prefix_parse(
                        first,
                        boost::parser::detail::text::null_sentinel,
                        parser);
                BOOST_TEST(chars);
                BOOST_TEST(*chars == std::vector<std::string>({"yay"}));
            }
        }
        {
            std::string str = "yay,";
            {
                std::vector<std::string> chars;
                BOOST_TEST(!parse(str, parser, chars));
            }
            {
                std::vector<std::string> chars;
                auto first = str.c_str();
                BOOST_TEST(prefix_parse(
                    first,
                    boost::parser::detail::text::null_sentinel,
                    parser,
                    chars));
                BOOST_TEST(chars == std::vector<std::string>({"yay"}));
            }

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser);
                BOOST_TEST(!chars);
            }
            {
                auto first = str.c_str();
                std::optional<std::vector<std::string>> const chars =
                    prefix_parse(
                        first,
                        boost::parser::detail::text::null_sentinel,
                        parser);
                BOOST_TEST(chars);
                BOOST_TEST(*chars == std::vector<std::string>({"yay"}));
            }
        }
        {
            char const * str = "yay,yay,yay";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<std::string>({"yay", "yay", "yay"}));

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(null_term(str), parser);
                BOOST_TEST(chars);
                BOOST_TEST(
                    *chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }
        }
    }

    {
        constexpr auto parser = string("yay") % ',';
        {
            char const * str = "";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(null_term(str), parser, char_(' '), chars));
            BOOST_TEST(chars == std::vector<std::string>{});
        }

        {
            char const * str = "";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser, char_(' '));
            BOOST_TEST(!chars);
        }
        {
            char const * str = "z";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(null_term(str), parser, char_(' '), chars));
            BOOST_TEST(chars == std::vector<std::string>{});
        }
        {
            char const * str = "z";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser, char_(' '));
            BOOST_TEST(!chars);
        }
        {
            char const * str = ",";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(null_term(str), parser, char_(' '), chars));
            BOOST_TEST(chars == std::vector<std::string>{});
        }
        {
            char const * str = ",";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser, char_(' '));
            BOOST_TEST(!chars);
        }
        {
            char const * str = " ,yay";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(null_term(str), parser, char_(' '), chars));
            BOOST_TEST(chars == std::vector<std::string>{});
        }
        {
            char const * str = " ,yay";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser, char_(' '));
            BOOST_TEST(!chars);
        }
        {
            char const * str = ", yay";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(null_term(str), parser, char_(' '), chars));
            BOOST_TEST(chars == std::vector<std::string>{});
        }
        {
            char const * str = ", yay";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser, char_(' '));
            BOOST_TEST(!chars);
        }
        {
            char const * str = ",yay ";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(null_term(str), parser, char_(' '), chars));
            BOOST_TEST(chars == std::vector<std::string>{});
        }
        {
            char const * str = ",yay ";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser, char_(' '));
            BOOST_TEST(!chars);
        }

        {
            char const * str = " , yay ";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(null_term(str), parser, char_(' '), chars));
            BOOST_TEST(chars == std::vector<std::string>{});
        }
        {
            char const * str = " , yay ";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser, char_(' '));
            BOOST_TEST(!chars);
        }
        {
            char const * str = "yay";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, char_(' '), chars));
            BOOST_TEST(chars == std::vector<std::string>({"yay"}));
        }
        {
            char const * str = "yay";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser, char_(' '));
            BOOST_TEST(chars);
            BOOST_TEST(*chars == std::vector<std::string>({"yay"}));
        }
        {
            char const * str = "yayyay";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(null_term(str), parser, char_(' '), chars));
        }
        {
            std::string str = "yayyay";
            std::vector<std::string> chars;
            auto first = str.c_str();
            BOOST_TEST(prefix_parse(
                first,
                boost::parser::detail::text::null_sentinel,
                parser,
                char_(' '),
                chars));
            BOOST_TEST(chars == std::vector<std::string>({"yay"}));
        }
        {
            char const * str = "yayyay";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser, char_(' '));
            BOOST_TEST(!chars);
        }
        {
            std::string str = "yayyay";
            auto first = str.c_str();
            std::optional<std::vector<std::string>> const chars = prefix_parse(
                first,
                boost::parser::detail::text::null_sentinel,
                parser,
                char_(' '));
            BOOST_TEST(chars);
            BOOST_TEST(*chars == std::vector<std::string>({"yay"}));
        }
        {
            char const * str = "yay,";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(null_term(str), parser, char_(' '), chars));
        }
        {
            std::string str = "yay,";
            std::vector<std::string> chars;
            auto first = str.c_str();
            BOOST_TEST(prefix_parse(
                first,
                boost::parser::detail::text::null_sentinel,
                parser,
                char_(' '),
                chars));
            BOOST_TEST(chars == std::vector<std::string>({"yay"}));
        }
        {
            char const * str = "yay,";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser, char_(' '));
            BOOST_TEST(!chars);
        }
        {
            std::string str = "yay,";
            auto first = str.c_str();
            std::optional<std::vector<std::string>> const chars = prefix_parse(
                first,
                boost::parser::detail::text::null_sentinel,
                parser,
                char_(' '));
            BOOST_TEST(chars);
            BOOST_TEST(*chars == std::vector<std::string>({"yay"}));
        }
        {
            char const * str = "yay,yay,yay";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, char_(' '), chars));
            BOOST_TEST(chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
        {
            char const * str = "yay,yay,yay";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser, char_(' '));
            BOOST_TEST(chars);
            BOOST_TEST(*chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
        {
            char const * str = " yay,yay,yay";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, char_(' '), chars));
            BOOST_TEST(chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
        {
            char const * str = " yay,yay,yay";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser, char_(' '));
            BOOST_TEST(chars);
            BOOST_TEST(*chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
        {
            char const * str = "yay ,yay,yay";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, char_(' '), chars));
            BOOST_TEST(chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }

        {
            char const * str = "yay ,yay,yay";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser, char_(' '));
            BOOST_TEST(chars);
            BOOST_TEST(*chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
        {
            char const * str = "yay, yay,yay";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, char_(' '), chars));
            BOOST_TEST(chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
        {
            char const * str = "yay, yay,yay";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser, char_(' '));
            BOOST_TEST(chars);
            BOOST_TEST(*chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
        {
            char const * str = "yay,yay ,yay";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, char_(' '), chars));
            BOOST_TEST(chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }

        {
            char const * str = "yay,yay ,yay";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser, char_(' '));
            BOOST_TEST(chars);
            BOOST_TEST(*chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
        {
            char const * str = "yay,yay, yay";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, char_(' '), chars));
            BOOST_TEST(chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }

        {
            char const * str = "yay,yay, yay";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser, char_(' '));
            BOOST_TEST(chars);
            BOOST_TEST(*chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
        {
            char const * str = "yay,yay,yay ";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, char_(' '), chars));
            BOOST_TEST(chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }

        {
            char const * str = "yay,yay,yay ";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser, char_(' '));
            BOOST_TEST(chars);
            BOOST_TEST(*chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
        {
            char const * str = " yay , yay , yay ";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, char_(' '), chars));
            BOOST_TEST(chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }

        {
            char const * str = " yay , yay , yay ";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser, char_(' '));
            BOOST_TEST(chars);
            BOOST_TEST(*chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
        {
            char const * str = "yay, yay, yay";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, char_(' '), chars));
            BOOST_TEST(chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }

        {
            char const * str = "yay, yay, yay";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser, char_(' '));
            BOOST_TEST(chars);
            BOOST_TEST(*chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
    }
}

// lexeme
{
    {
        constexpr auto parser = lexeme[string("yay") % ','];

        {
            std::string str = "yay, yay, yay";
            {
                std::vector<std::string> chars;
                BOOST_TEST(!parse(str, parser, char_(' '), chars));
            }
            {
                std::vector<std::string> chars;
                auto first = str.c_str();
                BOOST_TEST(prefix_parse(
                    first,
                    boost::parser::detail::text::null_sentinel,
                    parser,
                    char_(' '),
                    chars));
                BOOST_TEST(chars == std::vector<std::string>({"yay"}));
            }

            {
                char const * str = "yay, yay, yay";
                std::optional<std::vector<std::string>> const chars =
                    parse(null_term(str), parser, char_(' '));
                BOOST_TEST(!chars);
            }
            {
                std::string str = "yay, yay, yay";
                auto first = str.c_str();
                std::optional<std::vector<std::string>> const chars =
                    prefix_parse(
                        first,
                        boost::parser::detail::text::null_sentinel,
                        parser,
                        char_(' '));
                BOOST_TEST(chars);
                BOOST_TEST(*chars == std::vector<std::string>({"yay"}));
            }
        }
        {
            std::string str = " yay, yay, yay";
            {
                std::vector<std::string> chars;
                BOOST_TEST(!parse(str, parser, char_(' '), chars));
            }
            {
                std::vector<std::string> chars;
                auto first = str.c_str();
                BOOST_TEST(prefix_parse(
                    first,
                    boost::parser::detail::text::null_sentinel,
                    parser,
                    char_(' '),
                    chars));
                BOOST_TEST(chars == std::vector<std::string>({"yay"}));
            }

            {
                char const * str = " yay, yay, yay";
                std::optional<std::vector<std::string>> const chars =
                    parse(null_term(str), parser, char_(' '));
                BOOST_TEST(!chars);
            }
            {
                std::string str = " yay, yay, yay";
                auto first = str.c_str();
                std::optional<std::vector<std::string>> const chars =
                    prefix_parse(
                        first,
                        boost::parser::detail::text::null_sentinel,
                        parser,
                        char_(' '));
                BOOST_TEST(chars);
                BOOST_TEST(*chars == std::vector<std::string>({"yay"}));
            }
        }
    }

    {
        constexpr auto parser = lexeme[skip[string("yay") % ',']];

        {
            char const * str = "yay, yay, yay";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, char_(' '), chars));
            BOOST_TEST(chars == std::vector<std::string>({"yay", "yay", "yay"}));

            {
                char const * str = "yay, yay, yay";
                std::optional<std::vector<std::string>> const chars =
                    parse(null_term(str), parser, char_(' '));
                BOOST_TEST(chars);
                BOOST_TEST(
                    *chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }
        }
        {
            char const * str = " yay, yay, yay";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, char_(' '), chars));
            BOOST_TEST(chars == std::vector<std::string>({"yay", "yay", "yay"}));

            {
                char const * str = " yay, yay, yay";
                std::optional<std::vector<std::string>> const chars =
                    parse(null_term(str), parser, char_(' '));
                BOOST_TEST(chars);
                BOOST_TEST(
                    *chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }
        }
    }
}

// skip
{
    {
        constexpr auto parser = skip(char_(' '))[string("yay") % ','];

        {
            char const * str = "yay, yay, yay";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
        {
            char const * str = "yay, yay, yay";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser);
            BOOST_TEST(chars);
            BOOST_TEST(*chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
        {
            char const * str = " yay, yay, yay";
            std::vector<std::string> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
        {
            char const * str = " yay, yay, yay";
            std::optional<std::vector<std::string>> const chars =
                parse(null_term(str), parser);
            BOOST_TEST(chars);
            BOOST_TEST(*chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
    }
}

// combined_seq_and_or
{
    {
        constexpr auto parser = char_('a') >> char_('b') >> char_('c') |
                                char_('x') >> char_('y') >> char_('z');
        using tup = tuple<char, char, char>;

        {
            char const * str = "abc";
            tuple<char, char, char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == tup('c', '\0', '\0'));
        }

        {
            char const * str = "abc";
            std::optional<std::string> const chars = parse(null_term(str), parser);
            BOOST_TEST(chars);
            BOOST_TEST(*chars == "abc");
        }

        {
            char const * str = "xyz";
            tuple<char, char, char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == tup('z', '\0', '\0'));
        }
    }

    {
        constexpr auto parser = char_('a') >> string("b") >> char_('c') |
                                char_('x') >> string("y") >> char_('z');
        {
            char const * str = "abc";
            std::string chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == "abc");
        }

        {
            char const * str = "abc";
            std::optional<std::string> const chars = parse(null_term(str), parser);
            BOOST_TEST(chars);
            BOOST_TEST(*chars == "abc");
        }

        {
            char const * str = "xyz";
            std::string chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == "xyz");
        }
    }

    {
        constexpr auto parser = char_('a') >> char_('b') >> char_('c') |
                                char_('x') >> char_('y') >> char_('z');
        using tup = tuple<char, char, char>;

        {
            char const * str = "abc";
            tuple<std::any, std::any, std::any> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
        }

        {
            char const * str = "xyz";
            tuple<char, char, char> chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == tup('z', '\0', '\0'));
        }
    }

    {
        constexpr auto parser = !char_('a');
        char const * str = "a";
        BOOST_TEST(!parse(null_term(str), parser));
    }

    {
        constexpr auto parser = &char_('a');
        char const * str = "a";
        BOOST_TEST(!parse(null_term(str), parser));
    }
    {
        constexpr auto parser = &char_('a');
        std::string str = "a";
        auto first = str.c_str();
        BOOST_TEST(prefix_parse(first, boost::parser::detail::text::null_sentinel, parser));
    }

    {
#if defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-shift-op-parentheses"
#endif
        constexpr auto parser = (char_('a') >> string("b") > char_('c')) |
                                (char_('x') >> string("y") >> char_('z'));
#if defined(__clang__)
#pragma GCC diagnostic pop
#endif
        {
            char const * str = "abc";
            std::string chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == "abc");
        }

        {
            char const * str = "abc";
            std::optional<std::string> const chars = parse(null_term(str), parser);
            BOOST_TEST(chars);
            BOOST_TEST(*chars == "abc");
        }

        {
            char const * str = "xyz";
            std::string chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == "xyz");
        }

        {
            char const * str = "abz";
            std::string chars;
            rethrow_error_handler eh;
            BOOST_TEST_THROWS(
                parse(null_term(str), with_error_handler(parser, eh), chars),
                parse_error<char const *>);
        }

        {
            char const * str = "abz";
            std::string chars;
            BOOST_TEST(!parse(null_term(str), parser, chars));
        }

        // stream_error_handler
        {
            char const * str = "ab";
            std::string chars;
            stream_error_handler eh;
            BOOST_TEST(!parse(null_term(str), with_error_handler(parser, eh), chars));
        }

        {
            char const * str = "abz";
            std::string chars;
            stream_error_handler eh("simple_parser.cpp");
            BOOST_TEST(!parse(null_term(str), with_error_handler(parser, eh), chars));
        }

        {
            char const * str = "ab";
            std::string chars;
            stream_error_handler eh("simple_parser.cpp");
            BOOST_TEST(!parse(null_term(str), with_error_handler(parser, eh), chars));
        }

        {
            char const * str = "abz";
            std::string chars;
            std::ostringstream oss;
            stream_error_handler eh("simple_parser.cpp", oss);
            BOOST_TEST(!parse(null_term(str), with_error_handler(parser, eh), chars));
            BOOST_TEST(
                oss.str() ==
                "simple_parser.cpp:1:2: error: Expected char_('c') "
                "here:\nabz\n  ^\n");
        }

        {
            char const * str = "abz";
            std::string chars;
            std::ostringstream err, warn;
            stream_error_handler eh("simple_parser.cpp", err, warn);
            BOOST_TEST(!parse(null_term(str), with_error_handler(parser, eh), chars));
            BOOST_TEST(
                err.str() ==
                "simple_parser.cpp:1:2: error: Expected char_('c') "
                "here:\nabz\n  ^\n");
            BOOST_TEST(warn.str() == "");
        }

        // callback_error_handler
        {
            {
                char const * str = "ab";
                std::string chars;
                callback_error_handler eh;
                BOOST_TEST(!parse(null_term(str), with_error_handler(parser, eh), chars));
            }

            std::string err_str;
            auto err = [&](std::string const & msg) { err_str = msg; };
            std::string warn_str;
            auto warn = [&](std::string const & msg) { warn_str = msg; };

            {
                char const * str = "abz";
                std::string chars;
                callback_error_handler eh(err);
                err_str = "";
                BOOST_TEST(!parse(null_term(str), with_error_handler(parser, eh), chars));
                BOOST_TEST(
                    err_str ==
                    "1:2: error: Expected char_('c') here:\nabz\n  ^\n");
            }

            {
                char const * str = "abz";
                std::string chars;
                callback_error_handler eh("simple_parser.cpp", err);
                err_str = "";
                BOOST_TEST(!parse(null_term(str), with_error_handler(parser, eh), chars));
                BOOST_TEST(
                    err_str ==
                    "simple_parser.cpp:1:2: error: Expected char_('c') "
                    "here:\nabz\n  ^\n");
            }

            {
                char const * str = "abz";
                std::string chars;
                callback_error_handler eh("simple_parser.cpp", err, warn);
                err_str = "";
                warn_str = "";
                BOOST_TEST(!parse(null_term(str), with_error_handler(parser, eh), chars));
                BOOST_TEST(
                    err_str ==
                    "simple_parser.cpp:1:2: error: Expected char_('c') "
                    "here:\nabz\n  ^\n");
                BOOST_TEST(warn_str == "");
            }
        }
    }

    {
#if defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-shift-op-parentheses"
#endif
        constexpr auto parser = (char_('a') >> string("b") > char_('c')) |
                                (char_('x') >> string("y") >> char_('z'));
#if defined(__clang__)
#pragma GCC diagnostic pop
#endif
        {
            char const * str = "abc";
            std::string chars;
            BOOST_TEST(parse(null_term(str), parser, chars));
            BOOST_TEST(chars == "abc");
        }

        {
            char const * str = "abc";
            std::optional<std::string> const chars = parse(null_term(str), parser);
            BOOST_TEST(chars);
            BOOST_TEST(*chars == "abc");
        }

        {
            char const * str = "xyz";
            std::string chars;
            BOOST_TEST(parse(null_term(str), parser, chars, trace::on));
            BOOST_TEST(chars == "xyz");
        }
    }
}

// broken_utf8
{
    constexpr char c = 0xcc;
    constexpr auto parser = *char_(c);
    {
        char const array[3] = {(char)0xcc, (char)0x80, 0}; // U+0300
        std::string str = array;
        std::string chars;
        auto first = str.begin();
        BOOST_TEST(prefix_parse(first, str.end(), parser, chars));
        char const expected[2] = {(char)0xcc, 0};
        BOOST_TEST(chars == expected); // Finds one match of the *char* 0xcc.
    }
#if defined(__cpp_char8_t)
    {
        std::u8string str = u8"\xcc\x80"; // U+0300
        std::string chars;
        auto first = str.begin();
        BOOST_TEST(prefix_parse(first, str.end(), parser, chars));
        BOOST_TEST(chars == ""); // Finds zero matches of the *code point* 0xcc.
    }
#endif
}

// attr_out_param_compat
{
    {
        namespace bp = boost::parser;
        auto const p = bp::string("foo");

        std::vector<char> result;
        bool const success = bp::parse("foo", p, result);
        BOOST_TEST(success && result == std::vector<char>({'f', 'o', 'o'}));
    }
    {
        namespace bp = boost::parser;
        auto const p = bp::string("foo");

        std::vector<int> result;
        bool const success = bp::parse("foo", p, result);
        BOOST_TEST(success && result == std::vector<int>({'f', 'o', 'o'}));
    }
    {
        namespace bp = boost::parser;
        auto const p = +(bp::cp - ' ') >> ' ' >> string("foo");

        using attr_type = decltype(bp::parse(u8"", p));
        static_assert(std::is_same_v<
                      attr_type,
                      std::optional<bp::tuple<std::string, std::string>>>);

        bp::tuple<std::vector<int>, std::string> result;
        bool const success = bp::parse(u8"rôle foo" | bp::as_utf8, p, result);
        using namespace bp::literals;

        assert(success);
        (void)success;
        assert(bp::get(result, 0_c) == std::vector<int>({'r', U'ô', 'l', 'e'}));
        assert(bp::get(result, 1_c) == "foo");
    }
    {
        namespace bp = boost::parser;
        auto const p = +(bp::cp - ' ') >> ' ' >> string("foo");

        using attr_type = decltype(bp::parse(u8"", p));
        static_assert(std::is_same_v<
                      attr_type,
                      std::optional<bp::tuple<std::string, std::string>>>);

        bp::tuple<std::vector<char>, std::string> result;
        bool const success = bp::parse(u8"rôle foo" | bp::as_utf8, p, result);
        using namespace bp::literals;

        assert(success);
        (void)success;
        // The 4 code points "rôle" get transcoded to 5 UTF-8 code points to fit in the std::string.
        assert(bp::get(result, 0_c) == std::vector<char>({'r', (char)0xc3, (char)0xb4, 'l', 'e'}));
        assert(bp::get(result, 1_c) == "foo");
    }
}

return boost::report_errors();
}
