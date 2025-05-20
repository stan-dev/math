// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/chrono.hpp>
#include <iostream>

int main()
{
    std::cout << boost::chrono::system_clock::now() << std::endl;
}
