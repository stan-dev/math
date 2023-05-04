// Copyright 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/system/result.hpp>
#include <boost/core/lightweight_test_trait.hpp>

using namespace boost::system;

struct X {};

int main()
{
    BOOST_TEST_TRAIT_SAME( result<int>::value_type, int );
    BOOST_TEST_TRAIT_SAME( result<X>::value_type, X );
    BOOST_TEST_TRAIT_SAME( result<void>::value_type, void );

    BOOST_TEST_TRAIT_SAME( result<int>::error_type, error_code );
    BOOST_TEST_TRAIT_SAME( result<int, X>::error_type, X );

    return boost::report_errors();
}
