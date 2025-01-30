//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale.hpp>
#include <iostream>
#include <set>
#include <string>
#include <vector>

using namespace boost::locale;

int main(int argc, char** argv)
{
    if(argc != 3) {
        std::cerr << "Usage backend locale\n";
        return 1;
    } else {
        boost::locale::localization_backend_manager mgr = boost::locale::localization_backend_manager::global();
        mgr.select(argv[1]);
        generator gen(mgr);
        // Set global locale to requested
        std::locale::global(gen(argv[2]));
    }

    // Read all strings into a vector
    std::vector<std::string> all;
    while(!std::cin.eof()) {
        std::string tmp;
        std::getline(std::cin, tmp);
        all.push_back(tmp);
    }
    for(int i = 0; i < 10000; i++) {
        std::vector<std::string> tmp = all;
        // std::locale can be used as object for comparison
        std::sort(tmp.begin(), tmp.end(), std::locale());
        if(i == 0) {
            for(const auto& s : tmp)
                std::cout << s << std::endl;
        }
    }
}
