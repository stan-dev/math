//
// Copyright (c) 2022 Alexander Grund
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale/config.hpp>

namespace boost { namespace locale { namespace test {
    // POSIX setenv/unsetenv for all platforms and with ISO C++ compiler setting
    int setenv(const char* key, const char* value);
    int unsetenv(const char* key);
}}} // namespace boost::locale::test
