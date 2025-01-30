// Copyright 2019 Peter Dimov
// Copyright 2022-2024 Antony Polukhin
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/stacktrace/stacktrace.hpp>
#include <iostream>

int main()
{
    std::cout << boost::stacktrace::stacktrace() << std::endl;
}
