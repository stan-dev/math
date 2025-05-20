// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/optional/optional.hpp>

int main()
{
    boost::optional<int> o1;
    boost::optional<int> o2( 0 );
    o1 = o2;
    return *o1;
}
