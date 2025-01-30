/*
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * https://www.boost.org/LICENSE_1_0.txt)
 *
 * Copyright (c) 2023 Andrey Semashev
 */
/*!
 * \file   self_contained_header.cpp
 * \author Andrey Semashev
 *
 * \brief  This file contains a test boilerplate for checking that every public header is self-contained and does not have any missing #includes.
 */

#define BOOST_SCOPE_TEST_INCLUDE_HEADER() <boost/scope/BOOST_SCOPE_TEST_HEADER>

#include BOOST_SCOPE_TEST_INCLUDE_HEADER()

int main(int, char*[])
{
    return 0;
}
