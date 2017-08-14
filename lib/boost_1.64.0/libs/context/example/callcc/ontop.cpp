
//          Copyright Oliver Kowalke 2016.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cstdlib>
#include <iostream>
#include <tuple>

#include <boost/context/continuation.hpp>

namespace ctx = boost::context;

int main() {
    ctx::continuation c;
    int data = 0;
    c = ctx::callcc( [](ctx::continuation && c) {
                        int data = c.get_data< int >();
                        std::cout << "f1: entered first time: " << data  << std::endl;
                        c = c.resume( data + 1);
                        data = c.get_data< int >();
                        std::cout << "f1: entered second time: " << data  << std::endl;
                        c = c.resume( data + 1);
                        data = c.get_data< int >();
                        std::cout << "f1: entered third time: " << data << std::endl;
                        return std::move( c);
                    },
                    data + 1);
    data = c.get_data< int >();
    std::cout << "f1: returned first time: " << data << std::endl;
    c = c.resume( data + 1);
    data = c.get_data< int >();
    std::cout << "f1: returned second time: " << data << std::endl;
    c = c.resume_with(
           [](ctx::continuation && c){
              int data = c.get_data< int >();
              std::cout << "f2: entered: " << data << std::endl;
              return data;
           },
           -1);
    std::cout << "f1: returned third time" << std::endl;
    std::cout << "main: done" << std::endl;
    return EXIT_SUCCESS;
}
