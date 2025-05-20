// Copyright (C) 2024 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//[ parsing_into_a_class_example
#include <boost/parser/parser.hpp>

#include <iostream>
#include <string>


namespace bp = boost::parser;

int main()
{
    std::cout << "Enter a string followed by two unsigned integers. ";
    std::string input;
    std::getline(std::cin, input);

    //[ parsing_into_a_class_str
    constexpr auto string_uint_uint =
        bp::lexeme[+(bp::char_ - ' ')] >> bp::uint_ >> bp::uint_;
    std::string string_from_parse;
    if (parse(input, string_uint_uint, bp::ws, string_from_parse))
        std::cout << "That yields this string: " << string_from_parse << "\n";
    else
        std::cout << "Parse failure.\n";
    //]

    std::cout << "Enter an unsigned integer followed by a string. ";
    std::getline(std::cin, input);
    std::cout << input << "\n";

    //[ parsing_into_a_class_vec_of_strs
    constexpr auto uint_string = bp::uint_ >> +bp::char_;
    std::vector<std::string> vector_from_parse;
    if (parse(input, uint_string, bp::ws, vector_from_parse)) {
        std::cout << "That yields this vector of strings:\n";
        for (auto && str : vector_from_parse) {
            std::cout << "  '" << str << "'\n";
        }
    } else {
        std::cout << "Parse failure.\n";
    }
    //]
}
//]
