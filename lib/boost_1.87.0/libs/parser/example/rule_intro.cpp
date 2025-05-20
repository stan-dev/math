// Copyright (C) 2020 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//[ rule_intro_example
#include <boost/parser/parser.hpp>

#include <deque>
#include <iostream>
#include <string>


namespace bp = boost::parser;


//[ rule_intro_rule_definition
//[ rule_intro_rule_definition_rule
bp::rule<struct doubles, std::vector<double>> doubles = "doubles";
//]
//[ rule_intro_rule_definition_rule_def
auto const doubles_def = bp::double_ % ',';
//]
//[ rule_intro_rule_definition_macro
BOOST_PARSER_DEFINE_RULES(doubles);
//]
//]

int main()
{
    std::cout << "Please enter a list of doubles, separated by commas. ";
    std::string input;
    std::getline(std::cin, input);

    //[ rule_intro_parse_call
    auto const result = bp::parse(input, doubles, bp::ws);
    //]

    if (result) {
        std::cout << "You entered:\n";
        for (double x : *result) {
            std::cout << x << "\n";
        }
    } else {
        std::cout << "Parse failure.\n";
    }
}
//]
