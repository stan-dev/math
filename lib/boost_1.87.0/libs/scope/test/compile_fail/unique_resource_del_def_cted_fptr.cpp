/*
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * https://www.boost.org/LICENSE_1_0.txt)
 *
 * Copyright (c) 2023 Andrey Semashev
 */
/*!
 * \file   unique_resource_del_def_cted_fptr.cpp
 * \author Andrey Semashev
 *
 * \brief  This file tests that \c unique_resource with a function pointer
 *         deleter cannot be default-constructed.
 */

#include <boost/scope/unique_resource.hpp>

int main()
{
    boost::scope::unique_resource< int, void (*)(int) > ur;

    return 0;
}
