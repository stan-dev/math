// Copyright (C) 2018 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include <boost/parser/parser.hpp>

#include <boost/core/lightweight_test.hpp>


using namespace boost::parser;

constexpr rule<struct abc_def_tag, std::string> abc_def = "abc or def";
constexpr auto abc_def_def = string("abc") | string("def");
BOOST_PARSER_DEFINE_RULES(abc_def);

auto const fail = [](auto & ctx) { _pass(ctx) = false; };
constexpr rule<struct fail_abc_pass_def_tag, std::string> fail_abc_pass_def =
    "abc";
constexpr auto fail_abc_pass_def_def = string("abc")[fail] | string("def");
BOOST_PARSER_DEFINE_RULES(fail_abc_pass_def);

auto const attr_to_val = [](auto & ctx) { _val(ctx) = _attr(ctx); };
constexpr rule<struct action_copy_abc_def_tag, std::string>
    action_copy_abc_def = "abc or def";
constexpr auto action_copy_abc_def_def =
    string("abc")[attr_to_val] | string("def")[attr_to_val];
BOOST_PARSER_DEFINE_RULES(action_copy_abc_def);

auto const abc_value = [](auto & ctx) { _val(ctx) = "abc"; };
auto const def_value = [](auto & ctx) { _val(ctx) = "def"; };
constexpr rule<struct rev_abc_def_tag, std::string> rev_abc_def = "abc or def";
constexpr auto rev_abc_def_def =
    string("abc")[def_value] | string("def")[abc_value];
BOOST_PARSER_DEFINE_RULES(rev_abc_def);

auto const append_attr = [](auto & ctx) { _locals(ctx) += _attr(ctx); };
auto const locals_to_val = [](auto & ctx) {
    _val(ctx) = std::move(_locals(ctx));
};
rule<struct locals_abc_def_tag, std::string, std::string> const locals_abc_def =
    "abc or def";
auto locals_abc_def_def = -string("abc")[append_attr] >>
                          -string("def")[append_attr] >> eps[locals_to_val];
BOOST_PARSER_DEFINE_RULES(locals_abc_def);

// for alternate invocables

auto drop_result = [](auto & ctx) {};
auto auto_assign = [](auto & ctx) { return _attr(ctx); };
auto auto_assign_multi_string_1 = [](std::string const & str1,
                                     std::string const & str2) {
    return str1 + str2;
};
auto auto_assign_multi_string_2 = [](std::string && str1, std::string && str2) {
    return str1 + str2;
};

constexpr rule<struct str_rule_1_tag, std::string> str_rule_1 = "str_rule_1";
constexpr auto str_rule_1_def = +char_;
BOOST_PARSER_DEFINE_RULES(str_rule_1);

constexpr rule<struct str_rule_2_tag, std::string> str_rule_2 = "str_rule_2";
constexpr auto str_rule_2_def = (+char_)[drop_result];
BOOST_PARSER_DEFINE_RULES(str_rule_2);

constexpr rule<struct str_rule_3_tag, std::string> str_rule_3 = "str_rule_3";
constexpr auto str_rule_3_def = (+char_)[auto_assign];
BOOST_PARSER_DEFINE_RULES(str_rule_3);

constexpr rule<struct str_rule_6_tag, std::string> str_rule_6 = "str_rule_6";
constexpr auto str_rule_6_def =
    (+(char_ - ' ') >> ' ' >> +char_)[auto_assign_multi_string_1];
BOOST_PARSER_DEFINE_RULES(str_rule_6);

constexpr rule<struct str_rule_7_tag, std::string> str_rule_7 = "str_rule_7";
constexpr auto str_rule_7_def =
    (+(char_ - ' ') >> ' ' >> +char_)[auto_assign_multi_string_2];
BOOST_PARSER_DEFINE_RULES(str_rule_7);

int main()
{

// side_effects
{
    int i = 0;
    auto increment_i = [&i](auto & ctx) { ++i; };

    using no_attribute_return = decltype(parse("xyz", char_('a')[increment_i]));
    static_assert(std::is_same_v<no_attribute_return, bool>);

    {
        std::string const str = "xyz";
        BOOST_TEST(!parse(str, char_('a')[increment_i]));
        BOOST_TEST(i == 0);
        BOOST_TEST(!parse(str, char_('x')[increment_i]));
        BOOST_TEST(i == 1);
        auto first = str.c_str();
        BOOST_TEST(prefix_parse(
            first, boost::parser::detail::text::null_sentinel, char_('x')[increment_i]));
        BOOST_TEST(i == 2);
        BOOST_TEST(!parse(str, char_('a')[increment_i]));
        BOOST_TEST(i == 2);
        BOOST_TEST(!parse(str, char_('x')[increment_i]));
        BOOST_TEST(i == 3);
        first = str.c_str();
        BOOST_TEST(prefix_parse(
            first, boost::parser::detail::text::null_sentinel, char_('x')[increment_i]));
        BOOST_TEST(i == 4);
    }
}

// pass
{
    {
        std::string const str = "xyz";
        auto const result = parse(str, abc_def);
        BOOST_TEST(!result);
    }
    {
        std::string const str = "abc";
        auto const result = parse(str, abc_def);
        BOOST_TEST(result);
        BOOST_TEST(*result == "abc");
    }
    {
        std::string const str = "def";
        auto const result = parse(str, abc_def);
        BOOST_TEST(result);
        BOOST_TEST(*result == "def");
    }
    {
        std::string const str = "abc";
        auto const result = parse(str, fail_abc_pass_def);
        BOOST_TEST(!result);
    }
    {
        std::string const str = "def";
        auto const result = parse(str, fail_abc_pass_def);
        BOOST_TEST(result);
        BOOST_TEST(*result == "def");
    }
}

// val_attr
{
    {
        std::string const str = "abc";
        auto const result = parse(str, action_copy_abc_def);
        BOOST_TEST(result);
        BOOST_TEST(*result == "abc");
    }
    {
        std::string const str = "def";
        auto const result = parse(str, action_copy_abc_def);
        BOOST_TEST(result);
        BOOST_TEST(*result == "def");
    }
    {
        std::string const str = "abc";
        auto const result = parse(str, rev_abc_def);
        BOOST_TEST(result);
        BOOST_TEST(*result == "def");
    }
    {
        std::string const str = "def";
        auto const result = parse(str, rev_abc_def);
        BOOST_TEST(result);
        BOOST_TEST(*result == "abc");
    }
}

// locals
{
    {
        std::string const str = "";
        auto const result = parse(str, locals_abc_def);
        BOOST_TEST(result);
        BOOST_TEST(*result == "");
    }
    {
        std::string const str = "abc";
        auto const result = parse(str, locals_abc_def);
        BOOST_TEST(result);
        BOOST_TEST(*result == "abc");
    }
    {
        std::string const str = "abcdef";
        auto const result = parse(str, locals_abc_def);
        BOOST_TEST(result);
        BOOST_TEST(*result == "abcdef");
    }
    {
        std::string const str = "def";
        auto const result = parse(str, locals_abc_def);
        BOOST_TEST(result);
        BOOST_TEST(*result == "def");
    }
}

// alternate_invocables
{
    {
        auto result_1 = parse("some text", str_rule_1);
        BOOST_TEST(result_1);
        BOOST_TEST(*result_1 == "some text");

        auto result_2 = parse("some text", str_rule_2);
        BOOST_TEST(result_2);
        BOOST_TEST(*result_2 == "");

        auto result_3 = parse("some text", str_rule_3);
        BOOST_TEST(result_3);
        BOOST_TEST(*result_3 == "some text");

        auto result_6 = parse("some text", str_rule_6);
        BOOST_TEST(result_6);
        BOOST_TEST(*result_6 == "sometext");

        auto result_7 = parse("some text", str_rule_7);
        BOOST_TEST(result_7);
        BOOST_TEST(*result_7 == "sometext");
    }
}

return boost::report_errors();
}
