/*  Copyright 2020 Rene Rivera
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or http://www.boost.org/LICENSE_1_0.txt)
 */

/*
This program is a compile test for support of C++11. If it compiles
successfully some key parts of C++11 the B2 engine requires are
available. This is used by the build script to guess and check the
compiler to build the engine with.
*/

// Some headers we depend on..
#include <thread>


int main()
{
    // Check for basic thread calls.
    { auto _ = std::thread::hardware_concurrency(); }
}
