// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/smart_ptr/detail/sp_type_traits.hpp>
#include <boost/core/lightweight_test_trait.hpp>

struct X;

int main()
{
    using boost::detail::sp_is_bounded_array;

    BOOST_TEST_TRAIT_FALSE(( sp_is_bounded_array<void> ));
    BOOST_TEST_TRAIT_FALSE(( sp_is_bounded_array<int> ));
    BOOST_TEST_TRAIT_FALSE(( sp_is_bounded_array<X> ));

    BOOST_TEST_TRAIT_FALSE(( sp_is_bounded_array<int[]> ));
    BOOST_TEST_TRAIT_FALSE(( sp_is_bounded_array<X[]> ));

    BOOST_TEST_TRAIT_TRUE(( sp_is_bounded_array<int[1]> ));
    BOOST_TEST_TRAIT_TRUE(( sp_is_bounded_array<int[7]> ));
    BOOST_TEST_TRAIT_TRUE(( sp_is_bounded_array<X[1]> ));
    BOOST_TEST_TRAIT_TRUE(( sp_is_bounded_array<X[7]> ));

    return boost::report_errors();
}
