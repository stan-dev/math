
// Copyright 2017 Daniel James.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Not a test, just a small program to write out configuration info

#include <boost/container_hash/hash.hpp>
#include <iostream>
#include <algorithm>
#include <limits>
#include <climits>

#if defined(BOOST_MSVC)

struct msvc_version {
    unsigned version;
    char const* description;

    friend bool operator<(msvc_version const& v1, msvc_version const& v2) {
        return v1.version < v2.version;
    }
};

void write_compiler_info() {
    // From:
    // https://en.wikipedia.org/wiki/Microsoft_Visual_C%2B%2B
    // https://blogs.msdn.microsoft.com/vcblog/2017/11/15/side-by-side-minor-version-msvc-toolsets-in-visual-studio-2017/
    msvc_version versions[] = {
        {0, "Old Visual C++"},
        {1000, "Visual C++ 4.x, VS4.0?"},
        {1100, "Visual C++ 5.0, VS97"},
        {1200, "Visual C++ 6.0, VS6.0"},
        {1300, "Visual C++ 7.0, VS.NET 2002"},
        {1310, "Visual C++ 7.1, VS.NET 2003"},
        {1400, "Visual C++ 8.0, VS2005"},
        {1500, "Visual C++ 9.0, VS2008"},
        {1600, "Visual C++ 10.0, VS2010"},
        {1700, "Visual C++ 11.0, VS2012"},
        {1800, "Visual C++ 12.0, VS2013"},
        {1900, "Visual C++ 14.00, VS2015"},
        {1910, "Visual C++ 14.1x, VS2017"},
        {1920, "Visual C++ 14.2x, VS2019"},
        {1930, "Visual C++ 14.3x, VS2022"},
    };

    msvc_version msvc = { BOOST_MSVC, "" };
    msvc_version* v = std::upper_bound(versions,
        versions + sizeof(versions) / sizeof(*versions),
        msvc) - 1;
    unsigned difference = msvc.version - v->version;

    std::cout << v->description;
    if (difference) {
        std::cout << " +" << difference;
    }
    std::cout << std::endl;
}

#else

void write_compiler_info() {
}

#endif

#define PRINT(x) std::cout << #x ": " << x << std::endl

int main() {
    write_compiler_info();
    std::cout << std::endl;

    PRINT(__cplusplus);
    PRINT(BOOST_CXX_VERSION);
    std::cout << std::endl;

#if defined(BOOST_NO_CXX17_HDR_STRING_VIEW)
    std::cout << "No <string_view>" << std::endl;
#else
    std::cout << "Has <string_view>" << std::endl;
#endif

#if defined(BOOST_NO_CXX17_HDR_OPTIONAL)
    std::cout << "No <optional>" << std::endl;
#else
    std::cout << "Has <optional>" << std::endl;
#endif

#if defined(BOOST_NO_CXX17_HDR_VARIANT)
    std::cout << "No <variant>" << std::endl;
#else
    std::cout << "Has <variant>" << std::endl;
#endif

#if defined(BOOST_NO_CXX11_HDR_TYPEINDEX)
    std::cout << "No <typeindex>" << std::endl;
#else
    std::cout << "Has <typeindex>" << std::endl;
#endif

#if defined(BOOST_NO_CXX11_HDR_SYSTEM_ERROR)
    std::cout << "No <system_error>" << std::endl;
#else
    std::cout << "Has <system_error>" << std::endl;
#endif

    std::cout << std::endl;

    PRINT(CHAR_BIT);
    std::cout << std::endl;

    PRINT(sizeof(std::size_t)*CHAR_BIT);
    PRINT(std::numeric_limits<std::size_t>::digits);
    std::cout << std::endl;

    PRINT(sizeof(float)*CHAR_BIT);
    PRINT(std::numeric_limits<float>::digits);
    std::cout << std::endl;

    PRINT(sizeof(double)*CHAR_BIT);
    PRINT(std::numeric_limits<double>::digits);
    std::cout << std::endl;

    PRINT(sizeof(long double)*CHAR_BIT);
    PRINT(std::numeric_limits<long double>::digits);
}
