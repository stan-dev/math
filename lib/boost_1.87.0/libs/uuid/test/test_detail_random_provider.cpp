//
// Copyright (c) 2017 James E. King III
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
//   https://www.boost.org/LICENSE_1_0.txt)
//
// Positive and negative testing for detail::random_provider
//

#include <boost/uuid/detail/random_provider.hpp>
#include <boost/uuid/entropy_error.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/array.hpp>
#include <boost/core/lightweight_test.hpp>
#include <limits>
#include <cstdint>
#include <string.h>

int main()
{
    boost::uuids::detail::random_provider provider;
    boost::array<unsigned int, 2> ints;

    // test generator()
    ints[0] = 0;
    ints[1] = 0;
    provider.generate(ints.begin(), ints.end());
    BOOST_TEST_NE(ints[0], ints[1]);

    // test name()
    BOOST_TEST_GT(strlen(provider.name()), 4u);

    // test the equivalent of get_random_bytes()
    std::uint32_t buf1[16];
    std::uint32_t buf2[16];
    provider.generate(buf1, buf1 + 16);
    provider.generate(buf2, buf2 + 16);
    BOOST_TEST_NE(0, memcmp(buf1, buf2, 64));

    return boost::report_errors();
}
