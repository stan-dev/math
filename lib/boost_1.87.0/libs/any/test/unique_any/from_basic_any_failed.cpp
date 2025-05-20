//  Copyright Antony Polukhin, 2013-2024.
//
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt).

#include <boost/any/basic_any.hpp>
#include <boost/any/unique_any.hpp>

int main() {
    boost::anys::basic_any<> a;
    boost::anys::unique_any b(a);
    (void)b;
}
