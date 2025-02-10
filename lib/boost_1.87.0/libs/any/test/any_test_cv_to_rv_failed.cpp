//  Unit test for boost::any.
//
//  See http://www.boost.org for most recent version, including documentation.
//
//  Copyright Antony Polukhin, 2013-2024.
//
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt).

#include <cstdlib>
#include <string>
#include <utility>

#include <boost/any.hpp>

int main() {
    boost::any const cvalue(10);
    int i = boost::any_cast<int&&>(cvalue);
    (void)i;
}

