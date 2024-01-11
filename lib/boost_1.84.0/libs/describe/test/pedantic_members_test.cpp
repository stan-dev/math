// Copyright 2022 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/describe/members.hpp>
#include <boost/core/lightweight_test_trait.hpp>

struct X {};

int main()
{
#if defined(BOOST_DESCRIBE_CXX11)

    BOOST_TEST_TRAIT_FALSE((boost::describe::has_describe_members<X>));

#endif

    return boost::report_errors();
}
