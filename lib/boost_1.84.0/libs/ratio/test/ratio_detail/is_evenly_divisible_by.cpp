// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/ratio/detail/is_evenly_divisible_by.hpp>
#include <boost/ratio/ratio.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <cstdint>

using boost::ratio_detail::is_evenly_divisible_by;
using boost::ratio;

int main()
{
    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<1, 1>, ratio<1, 1> >));
    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<2, 1>, ratio<1, 1> >));
    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<2, 1>, ratio<2, 1> >));
    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<2, 1>, ratio<1, 2> >));

    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<INTMAX_MAX, 1>, ratio<1, 1> >));
    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<INTMAX_MAX, 1>, ratio<INTMAX_MAX, 1> >));
    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<INTMAX_MAX, 1>, ratio<1, INTMAX_MAX> >));

    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<4, 1>, ratio<2, 1> >));
    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<1, 2>, ratio<1, 4> >));

    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<15, 2>, ratio<3, 4> >));

    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<INTMAX_MAX / 2 * 2, 1>, ratio<INTMAX_MAX / 2, 1> >));
    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<1, INTMAX_MAX / 2>, ratio<1, INTMAX_MAX / 2 * 2> >));

    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<0, 1>, ratio<INTMAX_MAX, 1> >));

    BOOST_TEST_TRAIT_FALSE((is_evenly_divisible_by< ratio<1, 1>, ratio<2, 1> >));
    BOOST_TEST_TRAIT_FALSE((is_evenly_divisible_by< ratio<1, 1>, ratio<0, 1> >));

    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<-1, 1>, ratio<1, 1> >));
    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<1, 1>, ratio<-1, 1> >));
    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<-2, 1>, ratio<1, 1> >));
    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<2, 1>, ratio<-1, 1> >));
    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<-2, 1>, ratio<2, 1> >));
    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<2, 1>, ratio<-2, 1> >));
    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<-2, 1>, ratio<1, 2> >));
    BOOST_TEST_TRAIT_TRUE((is_evenly_divisible_by< ratio<2, 1>, ratio<-1, 2> >));

    return boost::report_errors();
}
