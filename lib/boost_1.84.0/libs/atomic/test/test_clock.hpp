//  Copyright (c) 2020 Andrey Semashev
//
//  Distributed under the Boost Software License, Version 1.0.
//  See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_ATOMIC_TEST_TEST_CLOCK_HPP_INCLUDED_
#define BOOST_ATOMIC_TEST_TEST_CLOCK_HPP_INCLUDED_

#include <boost/config.hpp>
#if defined(BOOST_WINDOWS)
#include <boost/winapi/config.hpp>
#include <boost/winapi/basic_types.hpp>
#include <boost/winapi/time.hpp>
#include <ratio>
#endif
#include <chrono>

namespace chrono = std::chrono;

#if defined(BOOST_LIBSTDCXX_VERSION) && BOOST_LIBSTDCXX_VERSION < 40700
typedef chrono::monotonic_clock steady_clock;
#else
typedef chrono::steady_clock steady_clock;
#endif

#if defined(BOOST_WINDOWS)

// On Windows high precision clocks tend to cause spurious test failures because threads wake up earlier than expected.
// Use a lower precision steady clock for tests.
struct test_clock
{
#if BOOST_USE_WINAPI_VERSION >= BOOST_WINAPI_VERSION_WIN6
    typedef boost::winapi::ULONGLONG_ rep;
#else
    typedef boost::winapi::DWORD_ rep;
#endif
    typedef std::milli period;
    typedef chrono::duration< rep, period > duration;
    typedef chrono::time_point< test_clock, duration > time_point;

    static BOOST_CONSTEXPR_OR_CONST bool is_steady = true;

    static time_point now() BOOST_NOEXCEPT
    {
#if BOOST_USE_WINAPI_VERSION >= BOOST_WINAPI_VERSION_WIN6
        rep ticks = boost::winapi::GetTickCount64();
#else
        rep ticks = boost::winapi::GetTickCount();
#endif
        return time_point(duration(ticks));
    }
};

#else
typedef steady_clock test_clock;
#endif

#endif // BOOST_ATOMIC_TEST_TEST_CLOCK_HPP_INCLUDED_
