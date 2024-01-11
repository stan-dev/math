// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/ratio/detail/is_ratio.hpp>
#include <boost/ratio/ratio.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <cstdint>

using boost::ratio_detail::is_ratio;
using boost::ratio;

int main()
{
    BOOST_TEST_TRAIT_TRUE((is_ratio< ratio<1, 1> >));
    BOOST_TEST_TRAIT_TRUE((is_ratio< ratio<INTMAX_MAX, 1> >));
    BOOST_TEST_TRAIT_TRUE((is_ratio< ratio<1, INTMAX_MAX> >));

    BOOST_TEST_TRAIT_FALSE((is_ratio< int >));
    BOOST_TEST_TRAIT_FALSE((is_ratio< void >));

    return boost::report_errors();
}
