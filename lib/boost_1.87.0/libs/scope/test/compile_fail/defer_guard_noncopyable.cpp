/*
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * https://www.boost.org/LICENSE_1_0.txt)
 *
 * Copyright (c) 2023 Andrey Semashev
 */
/*!
 * \file   defer_guard_noncopyable.cpp
 * \author Andrey Semashev
 *
 * \brief  This file tests that \c defer_guard is noncopyable.
 */

#include <boost/scope/defer.hpp>
#include "function_types.hpp"

int main()
{
    int n = 0;
    boost::scope::defer_guard< normal_func > guard1{ normal_func(n) };
    boost::scope::defer_guard< normal_func > guard2 = guard1;

    return 0;
}
