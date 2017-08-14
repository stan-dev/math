
//          Copyright Oliver Kowalke 2016.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cstdlib>
#include <iostream>

#include <boost/context/continuation.hpp>

namespace ctx = boost::context;

ctx::continuation f1( ctx::continuation && c) {
    int data = c.get_data< int >();
    std::cout << "f1: entered first time: " << data << std::endl;
    c = c.resume( data + 2);
    data = c.get_data< int >();
    std::cout << "f1: entered second time: " << data << std::endl;
    return std::move( c);
}

int main() {
    ctx::continuation c;
    int data = 1;
    c = ctx::callcc( f1, data);
    data = c.get_data< int >();
    std::cout << "f1: returned first time: " << data << std::endl;
    c = c.resume( data + 2);
    if ( c.data_available() ) {
        data = c.get_data< int >();
        std::cout << "f1: returned second time: " << data << std::endl;
    } else {
        std::cout << "f1: returned second time: no data" << std::endl;
    }
    std::cout << "main: done" << std::endl;
    return EXIT_SUCCESS;
}
