/*
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * https://www.boost.org/LICENSE_1_0.txt)
 *
 * Copyright (c) 2024 Andrey Semashev
 */
/*!
 * \file   unique_resource_bad_rtraits_mkdef_nonassignable.cpp
 * \author Andrey Semashev
 *
 * \brief  This file tests that \c unique_resource rejects resource
 *         traits where the resource cannot be assigned from the
 *         result of `make_default`.
 */

#include <boost/scope/unique_resource.hpp>

struct resource
{
    resource(int) noexcept {}
    resource& operator= (int) = delete;
};

struct resource_deleter
{
    using result_type = void;
    result_type operator() (resource const&) const noexcept {}
};

struct bad_resource_traits
{
    static int make_default() noexcept { return 10; }
    static bool is_allocated(resource const& res) noexcept { return false; }
};

int main()
{
    boost::scope::unique_resource< resource, resource_deleter, bad_resource_traits > ur;

    return 0;
}
