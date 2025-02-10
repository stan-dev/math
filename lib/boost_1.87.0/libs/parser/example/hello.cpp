// Copyright (C) 2020 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//[ hello_example
#include <boost/parser/parser.hpp>

#include <iostream>
#include <string>


namespace bp = boost::parser;

int main(int argc, char const * argv[])
{
    std::string input = "World";
    if (1 < argc)
        input = argv[1];

    std::string result;
    bp::parse(input, *bp::char_, result);
    std::cout << "Hello, " << result << "!\n";
}
//]
