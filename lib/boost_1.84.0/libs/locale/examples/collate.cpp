//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <iostream>
#include <set>
#include <string>

#include <boost/locale.hpp>

using namespace boost::locale;

int main()
{
    generator gen;
    // Set global locale
    std::locale::global(gen(""));

    // Create a set that includes all strings sorted in alphabetical order
    // std::locale can be used as object for comparison
    std::set<std::string, std::locale> all_strings;

    // Read all strings into the set
    while(!std::cin.eof()) {
        std::string tmp;
        std::getline(std::cin, tmp);
        all_strings.insert(tmp);
    }
    // Print them out
    for(const std::string& str : all_strings)
        std::cout << str << std::endl;
}
