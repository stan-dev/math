/*
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * https://www.boost.org/LICENSE_1_0.txt)
 *
 * Copyright (c) 2023 Andrey Semashev
 */
/*!
 * \file   scope_exit_cond_def_cted_fptr.cpp
 * \author Andrey Semashev
 *
 * \brief  This file tests that \c scope_exit with a function pointer
 *         condition function object cannot be default-constructed.
 */

#include <boost/scope/scope_exit.hpp>
#include "function_types.hpp"

int main()
{
    int n = 0;
    boost::scope::scope_exit< normal_func, bool (*)() > guard{ normal_func(n) };

    return 0;
}
