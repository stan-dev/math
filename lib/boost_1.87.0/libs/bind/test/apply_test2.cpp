// Copyright 2021 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/bind/apply.hpp>
#include <boost/core/lightweight_test_trait.hpp>

int main()
{
    BOOST_TEST_TRAIT_SAME(void, boost::apply<void>::result_type);
    BOOST_TEST_TRAIT_SAME(int&, boost::apply<int&>::result_type);

    return boost::report_errors();
}
