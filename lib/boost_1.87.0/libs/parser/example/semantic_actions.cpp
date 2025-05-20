// Copyright (C) 2020 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//[ semantic_action_example
#include <boost/parser/parser.hpp>

#include <iostream>
#include <string>


namespace bp = boost::parser;

int main()
{
    std::cout << "Enter a list of doubles, separated by commas. ";
    std::string input;
    std::getline(std::cin, input);

    std::vector<double> result;
//[ semantic_action_example_lambda
    auto const action = [&result](auto & ctx) {
        std::cout << "Got one!\n";
        result.push_back(_attr(ctx));
    };
//]
    auto const action_parser = bp::double_[action];
    auto const success = bp::parse(input, action_parser % ',', bp::ws);

    if (success) {
        std::cout << "You entered:\n";
        for (double x : result) {
            std::cout << x << "\n";
        }
    } else {
        std::cout << "Parse failure.\n";
    }
}
//]
