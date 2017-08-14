
//          Copyright Oliver Kowalke 2016.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cstdlib>
#include <iostream>
#include <string>

#include <boost/context/continuation.hpp>

namespace ctx = boost::context;

ctx::continuation f1( ctx::continuation && c) {
    int i = 3;
    c = c.resume( i);
    std::string s{ "abc" };
    c = c.resume( s);
    i = 7; s = "xyz";
    c = c.resume( i, s);
    c = c.resume();
    return std::move( c);
}

int main() {
    ctx::continuation c = ctx::callcc( f1);
    int i = c.get_data< int >();
    std::cout << "f1: returned : " << i << std::endl;
    c = c.resume();
    std::string s = c.get_data< std::string >();
    std::cout << "f1: returned : " << s << std::endl;
    c = c.resume();
    std::tie(i,s)=c.get_data< int, std::string >();
    std::cout << "f1: returned : " << i << ", " << s << std::endl;
    c = c.resume();
    std::cout << std::boolalpha;
    std::cout << "f1: returned data : " << c.data_available() << std::endl;
    std::cout << "main: done" << std::endl;
    return EXIT_SUCCESS;
}
