// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/cppalliance/charconv/issues/156

#include <boost/charconv.hpp>
#include <boost/core/lightweight_test.hpp>
#include <random>
#include <iostream>

#ifdef BOOST_CHARCONV_HAS_BRAINFLOAT16

int main()
{
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<float> dist (1e10F, 1e12F);

    for (int i = 0; i < 1024; ++i)
    {
        char buffer[256] {};
        auto val = static_cast<std::bfloat16_t>(dist(rng));
        auto r = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), val, boost::charconv::chars_format::scientific);
        if (!(BOOST_TEST(r) && BOOST_TEST(std::strlen(buffer) <= 9))) // 4 digits prec + period + e+ + 2 digits of exp
        {
            std::cerr << "Incorrect for: " << val << " returned: " << buffer << std::endl; // LCOV_EXCL_LINE
        }
    }

    return boost::report_errors();
}

#else

int main()
{
    return 0;
}

#endif
