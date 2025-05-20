// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/smart_ptr/detail/sp_type_traits.hpp>
#include <boost/core/lightweight_test_trait.hpp>

struct X;

int main()
{
    using boost::detail::sp_type_identity;

    BOOST_TEST_TRAIT_SAME( sp_type_identity<void>::type, void );
    BOOST_TEST_TRAIT_SAME( sp_type_identity<int>::type, int );
    BOOST_TEST_TRAIT_SAME( sp_type_identity<X>::type, X );

    BOOST_TEST_TRAIT_SAME( sp_type_identity<void const>::type, void const );
    BOOST_TEST_TRAIT_SAME( sp_type_identity<int const>::type, int const );
    BOOST_TEST_TRAIT_SAME( sp_type_identity<X const>::type, X const );

    BOOST_TEST_TRAIT_SAME( sp_type_identity<void(int, X)>::type, void(int, X) );

    BOOST_TEST_TRAIT_SAME( sp_type_identity<int[]>::type, int[] );
    BOOST_TEST_TRAIT_SAME( sp_type_identity<X[]>::type, X[] );

    BOOST_TEST_TRAIT_SAME( sp_type_identity<int[1]>::type, int[1] );
    BOOST_TEST_TRAIT_SAME( sp_type_identity<X[2]>::type, X[2] );

    return boost::report_errors();
}
