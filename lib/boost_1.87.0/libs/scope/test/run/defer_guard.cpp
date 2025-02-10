/*
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * https://www.boost.org/LICENSE_1_0.txt)
 *
 * Copyright (c) 2023-2024 Andrey Semashev
 */
/*!
 * \file   defer_guard.cpp
 * \author Andrey Semashev
 *
 * \brief  This file contains tests for \c defer_guard.
 */

#include <boost/scope/defer.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <boost/config.hpp>
#include <utility>
#include <stdexcept>
#include "function_types.hpp"

#if defined(_MSC_VER) && !defined(__clang__)
// warning C4702: unreachable code
// This warning is triggered by tests that unconditionally throw exception at some point
// and have code after that (e.g. parts of scope guard constructor and a check that verifies
// that the following code is not reached).
#pragma warning(disable: 4702)
#endif

int g_n = 0;

void check_normal()
{
    int n = 0;
    {
        boost::scope::defer_guard< normal_func > guard{ normal_func(n) };
    }
    BOOST_TEST_EQ(n, 1);

    n = 0;
    {
        normal_func func(n);
        boost::scope::defer_guard< normal_func& > guard(func);
    }
    BOOST_TEST_EQ(n, 1);

    struct local
    {
        static void raw_func()
        {
            ++g_n;
        }
    };

    g_n = 0;
    {
        boost::scope::defer_guard< void (&)() > guard(local::raw_func);
    }
    BOOST_TEST_EQ(g_n, 1);
}

void check_throw()
{
    int n = 0;
    try
    {
        boost::scope::defer_guard< normal_func > guard{ normal_func(n) };
        throw std::runtime_error("error");
    }
    catch (...) {}
    BOOST_TEST_EQ(n, 1);

    n = 0;
    try
    {
        boost::scope::defer_guard< throw_on_copy_func > guard{ throw_on_copy_func(n) };
        BOOST_ERROR("An exception is expected to be thrown by throw_on_copy_func");
    }
    catch (...) {}
    BOOST_TEST_EQ(n, 1);

    n = 0;
    try
    {
        boost::scope::defer_guard< throw_on_move_func > guard{ throw_on_move_func(n) };
    }
    catch (...)
    {
        BOOST_ERROR("An exception is not expected to be thrown by throw_on_move_func (copy ctor should be used)");
    }
    BOOST_TEST_EQ(n, 1);

    n = 0;
    bool scope_ended = false, exception_thrown = false, func_destroyed = false;
    try
    {
        boost::scope::defer_guard< throw_on_call_func > guard{ throw_on_call_func(n, func_destroyed) };
        func_destroyed = false;
        scope_ended = true;
    }
    catch (...)
    {
        exception_thrown = true;
    }
    BOOST_TEST_EQ(n, 1);
    BOOST_TEST(scope_ended);
    BOOST_TEST(exception_thrown);
    BOOST_TEST(func_destroyed);
}

void check_deduction()
{
#if !defined(BOOST_NO_CXX17_DEDUCTION_GUIDES)
    int n = 0;
    {
        boost::scope::defer_guard guard{ normal_func(n) };
        BOOST_TEST_TRAIT_SAME(decltype(guard), boost::scope::defer_guard< normal_func >);
    }
    BOOST_TEST_EQ(n, 1);

    n = 0;
    {
        boost::scope::defer_guard guard([&n] { ++n; });
    }
    BOOST_TEST_EQ(n, 1);

    struct local
    {
        static void raw_func()
        {
            ++g_n;
        }

#if !defined(BOOST_SCOPE_NO_CXX17_NOEXCEPT_FUNCTION_TYPES)
        static void raw_func_noexcept() noexcept
        {
            ++g_n;
        }
#endif
    };

    g_n = 0;
    {
        boost::scope::defer_guard guard{ local::raw_func };
        BOOST_TEST_TRAIT_SAME(decltype(guard), boost::scope::defer_guard< void (*)() >);
    }
    BOOST_TEST_EQ(g_n, 1);

#if !defined(BOOST_SCOPE_NO_CXX17_NOEXCEPT_FUNCTION_TYPES)
    g_n = 0;
    {
        boost::scope::defer_guard guard{ local::raw_func_noexcept };
        BOOST_TEST_TRAIT_SAME(decltype(guard), boost::scope::defer_guard< void (*)() noexcept >);
    }
    BOOST_TEST_EQ(g_n, 1);
#endif

    n = 0;
    {
        BOOST_SCOPE_DEFER [&n] { ++n; };
        BOOST_SCOPE_DEFER [&n] { ++n; };
    }
    BOOST_TEST_EQ(n, 2);
#endif // !defined(BOOST_NO_CXX17_DEDUCTION_GUIDES)
}

int main()
{
    check_normal();
    check_throw();
    check_deduction();

    return boost::report_errors();
}
