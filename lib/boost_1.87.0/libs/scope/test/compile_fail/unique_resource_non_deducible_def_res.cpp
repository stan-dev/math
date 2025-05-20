/*
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * https://www.boost.org/LICENSE_1_0.txt)
 *
 * Copyright (c) 2023 Andrey Semashev
 */
/*!
 * \file   unique_resource_non_deducible_def_res.cpp
 * \author Andrey Semashev
 *
 * \brief  This file tests that \c unique_resource template parameters are not deducible from \c default_resource keyword.
 */

#include <boost/scope/unique_resource.hpp>

struct nop_resource_deleter
{
    template< typename Resource >
    void operator() (Resource const& res) const noexcept
    {
    }
};

int main()
{
    boost::scope::unique_resource ur{ boost::scope::default_resource, nop_resource_deleter() };

    return 0;
}
