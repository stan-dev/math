//
// Copyright (c) 2022 Dmitry Arkhipov <grisumbras@yandex.ru)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/json
//

#define BOOST_JSON_SOURCE

#include <boost/json/detail/config.hpp>

#ifdef BOOST_WINDOWS
# include <windows.h>
#endif

// X11 macros
#define Opposite
#define Bool

#include <boost/json/src.hpp>

#include "test_suite.hpp"

namespace boost {
namespace json {

class intrusive_macros_test
{
public:
    void run()
    {
        // do nothing
    }
};

TEST_SUITE(intrusive_macros_test, "boost.json.intrusive_macros");

} // namespace json
} // namespace boost
