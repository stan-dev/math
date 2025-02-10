// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <quadmath.h>

int main()
{
    __float128 f = -2.0Q;
    f = fabsq(f);

    return 0;
}
