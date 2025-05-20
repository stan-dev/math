// Copyright 2020, 2022 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/lambda2.hpp>

int main()
{
    using namespace boost::lambda2;
    return (_1 + _2 * _3)( 1, 2, 3) == 1 + 2 * 3? 0: 1;
}
