// Copyright (C) 2018 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include <boost/parser/parser.hpp>

#include <boost/core/lightweight_test.hpp>


using namespace boost::parser;

rule<class symbol_rule, std::string_view> const symrule = "symbols";
symbols<std::string_view> rule_symbols;
auto const fwd_attr = [](auto & ctx) { _val(ctx) = _attr(ctx); };
auto symrule_def = rule_symbols[fwd_attr];
BOOST_PARSER_DEFINE_RULES(symrule);

int main()
{
// symbols_empty
{
    symbols<int> roman_numerals;
    symbols<std::string> named_strings;

    std::string const str = "I";
    BOOST_TEST(!parse(str, roman_numerals));
    BOOST_TEST(!parse(str, named_strings));
}

// symbols_simple
{
    symbols<int> const roman_numerals = {
        {"I", 1}, {"V", 5}, {"X", 10}, {"L", 50}, {"C", 100}};
    symbols<std::string> const named_strings = {
        {"I", "1"}, {"V", "5"}, {"X", "10"}, {"L", "50"}, {"C", "100"}};

    {
        auto const result = parse("I", roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 1);
    }
    {
        auto const result = parse("I", named_strings);
        BOOST_TEST(result);
        BOOST_TEST(*result == "1");
    }

    {
        auto const result = parse("L", roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 50);
    }
    {
        auto const result = parse("L", named_strings);
        BOOST_TEST(result);
        BOOST_TEST(*result == "50");
    }
}

// symbols_max_munch
{
    symbols<int> const roman_numerals = {
        {"I", 1},
        {"II", 2},
        {"III", 3},
        {"IV", 4},
        {"V", 5},
        {"VI", 6},
        {"VII", 7},
        {"VIII", 8},
        {"IX", 9},

        {"X", 10},
        {"XX", 20},
        {"XXX", 30},
        {"XL", 40},
        {"L", 50},
        {"LX", 60},
        {"LXX", 70},
        {"LXXX", 80},
        {"XC", 90},

        {"C", 100},
        {"CC", 200},
        {"CCC", 300},
        {"CD", 400},
        {"D", 500},
        {"DC", 600},
        {"DCC", 700},
        {"DCCC", 800},
        {"CM", 900},

        {"M", 1000}};

    {
        auto const result = parse("I", roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 1);
    }
    {
        auto const result = parse("II", roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 2);
    }
    {
        auto const result = parse("III", roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 3);
    }
    {
        auto const result = parse("IV", roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 4);
    }
    {
        auto const result = parse("V", roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 5);
    }
    {
        auto const result = parse("VI", roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 6);
    }
    {
        auto const result = parse("VII", roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 7);
    }
    {
        auto const result = parse("VIII", roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 8);
    }
    {
        auto const result = parse("IX", roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 9);
    }
}

// symbols_mutating
{
    symbols<int> roman_numerals;
    roman_numerals.insert_for_next_parse("I", 1);
    roman_numerals.insert_for_next_parse("V", 5);
    roman_numerals.insert_for_next_parse("X", 10);
    auto const add_numeral = [&roman_numerals](auto & context) {
        using namespace boost::parser::literals;
        char chars[2] = {get(_attr(context), 0_c), 0};
        roman_numerals.insert(context, chars, get(_attr(context), 1_c));
    };
    auto const numerals_parser = (char_ >> int_)[add_numeral] >> roman_numerals;

    {
        auto const result = parse("L50L", numerals_parser);
        BOOST_TEST(result);
        BOOST_TEST(*result == 50);
        BOOST_TEST(!parse("L", roman_numerals));
    }
    {
        auto const result = parse("C100C", numerals_parser);
        BOOST_TEST(result);
        BOOST_TEST(*result == 100);
        BOOST_TEST(!parse("C", roman_numerals));
    }
    {
        auto const result = parse("L50C", numerals_parser);
        BOOST_TEST(!result);
    }
}

// insert/erase/clear
{
    symbols<int> roman_numerals;
    roman_numerals.insert_for_next_parse("I", 1);
    roman_numerals.insert_for_next_parse("V", 5);
    roman_numerals.insert_for_next_parse("X", 10);

    auto const insert_numeral = [&roman_numerals](auto & context) {
        using namespace boost::parser::literals;
        char chars[2] = {get(_attr(context), 0_c), 0};
        roman_numerals.insert(context, chars, get(_attr(context), 1_c));
    };
    auto const erase_numeral = [&roman_numerals](auto & context) {
        using namespace boost::parser::literals;
        char chars[2] = {_attr(context), 0};
        roman_numerals.erase(context, chars);
    };
    auto const clear_numerals = [&roman_numerals](auto & context) {
        roman_numerals.clear(context);
    };

    auto const add_parser = ("add" >> char_ >> int_)[insert_numeral];
    auto const delete_parser = ("del" >> char_)[erase_numeral];
    auto const clear_parser = lit("clear")[clear_numerals];

    auto const next_insert_numeral = [&roman_numerals](auto & context) {
        using namespace boost::parser::literals;
        char chars[2] = {get(_attr(context), 0_c), 0};
        roman_numerals.insert_for_next_parse(
            context, chars, get(_attr(context), 1_c));
    };
    auto const next_erase_numeral = [&roman_numerals](auto & context) {
        using namespace boost::parser::literals;
        char chars[2] = {_attr(context), 0};
        roman_numerals.erase_for_next_parse(context, chars);
    };
    auto const next_clear_numerals = [&roman_numerals](auto & context) {
        roman_numerals.clear_for_next_parse(context);
    };

    auto const next_add_parser =
        ("next-add" >> char_ >> int_)[next_insert_numeral];
    auto const next_delete_parser = ("next-del" >> char_)[next_erase_numeral];
    auto const next_clear_parser = lit("next-clear")[next_clear_numerals];

    {
        // add only for this parse
        auto result = parse("addL50L", add_parser >> roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 50);

        result = parse("L", roman_numerals);
        BOOST_TEST(!result);
    }
    {
        // add only for the next parse
        auto result = parse("next-addL50L", next_add_parser >> roman_numerals);
        BOOST_TEST(!result);

        result = parse("L", roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 50);
    }
    {
        // add for both parses
        auto result = parse("addL50next-addL50L",
                            add_parser >> next_add_parser >> roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 50);

        result = parse("L", roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 50);
    }
    {
        // add for both parses
        auto result = parse("next-addL50addL50L",
                            next_add_parser >> add_parser >> roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 50);

        result = parse("L", roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 50);
    }

    {
        // delete nonexistent entry
        auto result = parse("delQ", delete_parser);
        BOOST_TEST(result);
    }
    {
        // delete nonexistent entry
        auto result = parse("next-delQ", next_delete_parser);
        BOOST_TEST(result);
    }

    {
        // delete only for this parse
        BOOST_TEST(*parse("V", roman_numerals) == 5);

        auto result = parse("delVV", delete_parser >> roman_numerals);
        BOOST_TEST(!result);

        result = parse("V", roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 5);
    }
    {
        // delete only for the next parse
        BOOST_TEST(*parse("V", roman_numerals) == 5);

        auto result = parse("next-delVV", next_delete_parser >> roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 5);

        result = parse("V", roman_numerals);
        BOOST_TEST(!result);
    }
    {
        // delete for both parses
        BOOST_TEST(*parse("X", roman_numerals) == 10);

        auto result = parse(
            "delXnext-delXX",
            delete_parser >> next_delete_parser >> roman_numerals);
        BOOST_TEST(!result);

        result = parse("X", roman_numerals);
        BOOST_TEST(!result);
    }
    {
        // add and then delete only for this parse
        auto result = parse(
            "addC100CdelCC",
            add_parser >> roman_numerals >> delete_parser >> roman_numerals);
        BOOST_TEST(!result);
    }
    {
        // add for this parse and then delete only for the next parse
        auto result = parse(
            "addC100Cnext-delCC",
            add_parser >> roman_numerals >> next_delete_parser >>
            roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == std::tuple(100, 100));
    }

    {
        // clear only for this parse
        BOOST_TEST(*parse("I", roman_numerals) == 1);
        BOOST_TEST(*parse("L", roman_numerals) == 50);

        BOOST_TEST(!parse("clearI", clear_parser >> roman_numerals));
        BOOST_TEST(!parse("clearL", clear_parser >> roman_numerals));

        BOOST_TEST(*parse("I", roman_numerals) == 1);
        BOOST_TEST(*parse("L", roman_numerals) == 50);
    }

    {
        // clear only for the next parse
        BOOST_TEST(*parse("I", roman_numerals) == 1);
        BOOST_TEST(*parse("L", roman_numerals) == 50);

        auto result = parse("next-clearI", next_clear_parser >> roman_numerals);
        BOOST_TEST(result);
        BOOST_TEST(*result == 1);

        BOOST_TEST(!parse("I", roman_numerals));
        BOOST_TEST(!parse("L", roman_numerals));
    }

    {
        // parse using symbols directly -- not using the table within a rule
        rule_symbols.clear_for_next_parse();
        rule_symbols.insert_for_next_parse("I", "one");
        rule_symbols.insert_for_next_parse("L", "50");

        auto result = parse("I", rule_symbols);
        BOOST_TEST(result);
        BOOST_TEST(*result == "one");

        result = parse("L", rule_symbols);
        BOOST_TEST(result);
        BOOST_TEST(*result == "50");

        BOOST_TEST(!parse("X", rule_symbols));
    }

    {
        // symbols within a rule
        rule_symbols.clear_for_next_parse();
        rule_symbols.insert_for_next_parse("foo", "foofie");
        rule_symbols.insert_for_next_parse("bar", "barrie");

        auto result = parse("foo", symrule);
        BOOST_TEST(result);
        BOOST_TEST(*result == "foofie");

        result = parse("bar", symrule);
        BOOST_TEST(result);
        BOOST_TEST(*result == "barrie");

        BOOST_TEST(!parse("X", symrule));
        BOOST_TEST(!parse("I", symrule));
    }

    {
        // symbols within a rule w/error handler
        rule_symbols.clear_for_next_parse();
        rule_symbols.insert_for_next_parse("foo", "foofie");
        rule_symbols.insert_for_next_parse("bar", "barrie");

        callback_error_handler error_handler(
            [](std::string_view m) { std::cout << m << "\n"; });
        auto parser = with_error_handler(symrule, error_handler);

        auto result = parse("foo", parser);
        BOOST_TEST(result);
        BOOST_TEST(*result == "foofie");

        result = parse("bar", parser);
        BOOST_TEST(result);
        BOOST_TEST(*result == "barrie");

        BOOST_TEST(!parse("baz", parser));
    }
}

return boost::report_errors();
}
