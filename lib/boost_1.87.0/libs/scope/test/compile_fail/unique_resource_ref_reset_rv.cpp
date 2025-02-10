/*
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * https://www.boost.org/LICENSE_1_0.txt)
 *
 * Copyright (c) 2023 Andrey Semashev
 */
/*!
 * \file   unique_resource_ref_reset_rv.cpp
 * \author Andrey Semashev
 *
 * \brief  This file tests that \c unique_resource for reference types does not accept rvalues in \c reset().
 */

#include <boost/scope/unique_resource.hpp>

template< typename Resource >
struct empty_resource_deleter
{
    void operator() (Resource const& res) const noexcept
    {
    }
};

int main()
{
    int res = 10;
    boost::scope::unique_resource< int const&, empty_resource_deleter< int const& > > ur1{ res, empty_resource_deleter< int const& >() };
    ur1.reset(20);

    return 0;
}
