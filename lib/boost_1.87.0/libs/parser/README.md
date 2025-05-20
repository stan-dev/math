# parser

This is a parser combinator library for C++.  As a quick example of use, here
is a complete program that parses one or more `double`s separated by commas,
ignoring whitespace:

```c++
#include <boost/parser/parser.hpp>

#include <iostream>
#include <string>


namespace bp = boost::parser;

int main()
{
    std::cout << "Enter a list of doubles, separated by commas.  No pressure. ";
    std::string input;
    std::getline(std::cin, input);

    auto const result = bp::parse(
        input, bp::double_ >> *(',' >> bp::double_), bp::ws);

    if (result) {
        std::cout << "Great! It looks like you entered:\n";
        for (double x : *result) {
            std::cout << x << "\n";
        }
    } else {
        std::cout
            << "Good job!  Please proceed to the recovery annex for cake.\n";
    }
}
```

This library is header-only, and has no Boost dependencies by default.

Features:

- Parsers that parse a variety of things.
- Combining operations that make complex parsers out of simpler ones.
- Multiple ways of getting data out of the parse, including via callbacks.
- Sentinel- and range-friendly.
- Very Unicode friendliness.
- Excellent error reporting, via diagnostics like those produced by GCC and Clang.
- Trace support for debugging your parsers.
- Clever hacks to make compile time errors easier to deal with.  (These are totally optional.)

This library first appeared in Boost 1.87.0

Master status:

[![Ubuntu](https://github.com/tzlaine/parser/actions/workflows/ubuntu.yml/badge.svg?branch=master)](https://github.com/tzlaine/parser/actions/workflows/ubuntu.yml)

[![Fedora](https://github.com/tzlaine/parser/actions/workflows/fedora.yml/badge.svg?branch=master)](https://github.com/tzlaine/parser/actions/workflows/fedora.yml)

[![Windows MSVC](https://github.com/tzlaine/parser/actions/workflows/windows.yml/badge.svg?branch=master)](https://github.com/tzlaine/parser/actions/workflows/windows.yml)

[![macos-12 - Clang 14](https://github.com/tzlaine/parser/actions/workflows/macos-12.yml/badge.svg?branch=master)](https://github.com/tzlaine/parser/actions/workflows/macos-12.yml)

Develop status:

[![Ubuntu](https://github.com/tzlaine/parser/actions/workflows/ubuntu.yml/badge.svg?branch=develop)](https://github.com/tzlaine/parser/actions/workflows/ubuntu.yml)

[![Fedora](https://github.com/tzlaine/parser/actions/workflows/fedora.yml/badge.svg?branch=develop)](https://github.com/tzlaine/parser/actions/workflows/fedora.yml)

[![Windows MSVC](https://github.com/tzlaine/parser/actions/workflows/windows.yml/badge.svg?branch=develop)](https://github.com/tzlaine/parser/actions/workflows/windows.yml)

[![macos-12 - Clang 14](https://github.com/tzlaine/parser/actions/workflows/macos-12.yml/badge.svg?branch=develop)](https://github.com/tzlaine/parser/actions/workflows/macos-12.yml)

[![License](https://img.shields.io/badge/license-boost-brightgreen.svg)](LICENSE_1_0.txt)
