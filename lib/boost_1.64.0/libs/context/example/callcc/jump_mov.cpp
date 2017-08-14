
//          Copyright Oliver Kowalke 2016.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cstdlib>
#include <iostream>

#include <boost/context/continuation.hpp>

namespace ctx = boost::context;

class moveable {
public:
    bool    state;
    int     value;

    moveable() :
        state( false),
        value( -1) {
        }

    moveable( int v) :
        state( true),                                                       
        value( v) {
        }

    moveable( moveable && other) :                       
        state( other.state),
        value( other.value) {  
            other.state = false;
            other.value = -1;
        }    

    moveable & operator=( moveable && other) {
        if ( this == & other) return * this;
        state = other.state;
        value = other.value;
        other.state = false;                                                                     
        other.value = -1;
        return * this;
    }

    moveable( moveable const& other) = delete;
    moveable & operator=( moveable const& other) = delete;
};

ctx::continuation f1( ctx::continuation && c) {
    moveable data = c.get_data< moveable >();
    c = c.resume( std::move( data) );
    return std::move( c);
}

int main() {
    ctx::continuation c;
    moveable data1{ 3 };
    c = ctx::callcc( std::allocator_arg, ctx::fixedsize_stack{}, f1, std::move( data1) );
    moveable data2;
    data2 = c.get_data< moveable >();
    std::cout << "main: done" << std::endl;
    return EXIT_SUCCESS;
}
