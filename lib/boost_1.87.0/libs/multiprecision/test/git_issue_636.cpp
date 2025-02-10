///////////////////////////////////////////////////////////////////////////////
//  Copyright 2024 Matt Borland. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  See: https://github.com/boostorg/multiprecision/issues/635

#include <boost/multiprecision/float128.hpp>
#include <type_traits>

int main()
{
    static_assert(std::is_trivially_copyable<boost::multiprecision::float128>::value, "Should be trivial now");
    return 0;
}
