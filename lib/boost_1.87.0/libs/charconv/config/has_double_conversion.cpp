// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <double-conversion/double-conversion.h>

int main()
{
    using namespace double_conversion;

    const int flags = 0;
    int processed;

    auto converter = StringToDoubleConverter(flags, 0.0, 0.0, "inf", "nan");

    auto result = converter.StringToDouble("0.0", 3, &processed);

    return static_cast<int>(result);
}
