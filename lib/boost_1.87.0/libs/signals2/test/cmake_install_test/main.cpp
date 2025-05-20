// Copyright 2024 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/signals2.hpp>

int f()
{
    return 4;
}

int main()
{
    boost::signals2::signal<int()> sig;

    sig.connect( f );

    return sig() == 4? 0: 1;
}
