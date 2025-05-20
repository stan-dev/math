// Copyright (C) 2020 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//[ roman_numeral_example
#include <boost/parser/parser.hpp>

#include <iostream>
#include <string>


namespace bp = boost::parser;

int main()
{
    std::cout << "Enter a number using Roman numerals. ";
    std::string input;
    std::getline(std::cin, input);

    //[ roman_numeral_symbol_tables
    bp::symbols<int> const ones = {
        {"I", 1},
        {"II", 2},
        {"III", 3},
        {"IV", 4},
        {"V", 5},
        {"VI", 6},
        {"VII", 7},
        {"VIII", 8},
        {"IX", 9}};

    bp::symbols<int> const tens = {
        {"X", 10},
        {"XX", 20},
        {"XXX", 30},
        {"XL", 40},
        {"L", 50},
        {"LX", 60},
        {"LXX", 70},
        {"LXXX", 80},
        {"XC", 90}};

    bp::symbols<int> const hundreds = {
        {"C", 100},
        {"CC", 200},
        {"CCC", 300},
        {"CD", 400},
        {"D", 500},
        {"DC", 600},
        {"DCC", 700},
        {"DCCC", 800},
        {"CM", 900}};
    //]

    //[ roman_numeral_actions
    int result = 0;
    auto const add_1000 = [&result](auto & ctx) { result += 1000; };
    auto const add = [&result](auto & ctx) { result += _attr(ctx); };
    //]

    //[ roman_numeral_parser
    using namespace bp::literals;
    auto const parser =
        *'M'_l[add_1000] >> -hundreds[add] >> -tens[add] >> -ones[add];
    //]

    if (bp::parse(input, parser) && result != 0)
        std::cout << "That's " << result << " in Arabic numerals.\n";
    else
        std::cout << "That's not a Roman number.\n";
}
//]
