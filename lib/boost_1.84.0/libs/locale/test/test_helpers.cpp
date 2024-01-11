//
// Copyright (c) 2022 Alexander Grund
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// setenv etc. are an extension, so we need this for Cygwin and MinGW
#if defined(__CYGWIN__) && !defined(_GNU_SOURCE)
// The setenv family of functions is an extension on Cygwin
#    define _GNU_SOURCE 1
#endif
#if defined(__MINGW32__) && defined(__STRICT_ANSI__)
// On MinGW-w64 the workaround is not needed and leads to warnings
#    include <_mingw.h>
#    ifndef __MINGW64_VERSION_MAJOR
#        undef __STRICT_ANSI__
#    endif
#endif

#include <boostLocale/test/test_helpers.hpp>
#include <cstdlib>
#ifdef BOOST_WINDOWS
#    include <list>
#    include <string>
#endif

namespace boost { namespace locale { namespace test {
#ifdef BOOST_WINDOWS
    // Needed as strings become part of the environment
    static std::list<std::string> env_values;

    int setenv(const char* key, const char* value)
    {
        env_values.push_back(key + std::string("=") + value);
        return _putenv(env_values.back().c_str());
    }

    int unsetenv(const char* key)
    {
        return setenv(key, "");
    }

#else
    int setenv(const char* key, const char* value)
    {
        return ::setenv(key, value, 1);
    }

    int unsetenv(const char* key)
    {
        return ::unsetenv(key);
    }
#endif
}}} // namespace boost::locale::test
