
// Copyright 2017 Peter Dimov.
//
// Distributed under the Boost Software License, Version 1.0.
//
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt

// See library home page at http://www.boost.org/libs/system

// Avoid spurious VC++ warnings
# define _CRT_SECURE_NO_WARNINGS

#include <boost/system/error_code.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/config.hpp>
#include <system_error>
#include <cerrno>
#include <string>
#include <cstdio>

static void test_generic_category()
{
    boost::system::error_category const & bt = boost::system::generic_category();
    std::error_category const & st = bt;

    BOOST_TEST_CSTR_EQ( bt.name(), st.name() );
}

static void test_system_category()
{
    boost::system::error_category const & bt = boost::system::system_category();
    std::error_category const & st = bt;

    BOOST_TEST_CSTR_EQ( bt.name(), st.name() );
}

int main()
{
    test_generic_category();
    test_system_category();

    return boost::report_errors();
}
