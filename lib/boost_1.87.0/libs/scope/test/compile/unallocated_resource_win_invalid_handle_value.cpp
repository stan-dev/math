/*
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * https://www.boost.org/LICENSE_1_0.txt)
 *
 * Copyright (c) 2024 Andrey Semashev
 */
/*!
 * \file   unallocated_resource_win_invalid_handle_value.cpp
 * \author Andrey Semashev
 *
 * \brief  This file tests that \c unallocated_resource supports
 *         \c INVALID_HANDLE_VALUE on Windows, as is described in the docs.
 *
 * This code relies on MSVC-specific behavior, which allows* \c INVALID_HANDLE_VALUE
 * to be specified in non-type template parameters.
 */

#include <boost/config.hpp>

#if defined(BOOST_WINDOWS) && defined(BOOST_MSVC) && \
    !defined(BOOST_NO_CXX17_FOLD_EXPRESSIONS) && !defined(BOOST_NO_CXX17_AUTO_NONTYPE_TEMPLATE_PARAMS)

#include <cstddef>
#include <windows.h>
#include <boost/core/functor.hpp>
#include <boost/scope/unique_resource.hpp>

int main()
{
    boost::scope::unique_resource<
        HANDLE,
        boost::core::functor< CloseHandle >,
        boost::scope::unallocated_resource< INVALID_HANDLE_VALUE, (HANDLE)NULL >
    > ur;
    (void)ur;

    return 0;
}

#else // defined(BOOST_WINDOWS) ...

#include <boost/config/pragma_message.hpp>

BOOST_PRAGMA_MESSAGE("Skipping test because it is not supported on the platform/compiler/C++ version.")

int main()
{
}

#endif // defined(BOOST_WINDOWS) ...
