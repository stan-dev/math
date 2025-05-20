//  Copyright Antony Polukhin, 2013-2024.
//
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt).

#include <boost/lexical_cast.hpp>
#include <string>

int main()
{
    volatile int i = 42;
    boost::lexical_cast<std::string>(i);
    return 1;
}

