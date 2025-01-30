// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/ratio/detail/gcd_lcm.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <cstdint>

using boost::ratio_detail::gcd;
using boost::ratio_detail::lcm;

int main()
{
    BOOST_TEST_EQ( (gcd<1, 1>::value), 1 );
    BOOST_TEST_EQ( (gcd<1, -1>::value), 1 );
    BOOST_TEST_EQ( (gcd<2, 2>::value), 2 );
    BOOST_TEST_EQ( (gcd<-2, -2>::value), 2 );
    BOOST_TEST_EQ( (gcd<INTMAX_MAX, INTMAX_MAX>::value), INTMAX_MAX );
    BOOST_TEST_EQ( (gcd<INTMAX_MAX, -INTMAX_MAX>::value), INTMAX_MAX );
    BOOST_TEST_EQ( (gcd<-INTMAX_MAX, -INTMAX_MAX>::value), INTMAX_MAX );
    BOOST_TEST_EQ( (gcd<7*11, 11*13>::value), 11 );

    BOOST_TEST_EQ( (lcm<1, 1>::value), 1 );
    BOOST_TEST_EQ( (lcm<1, -1>::value), 1 );
    BOOST_TEST_EQ( (lcm<2, 2>::value), 2 );
    BOOST_TEST_EQ( (lcm<-2, -2>::value), 2 );
    BOOST_TEST_EQ( (lcm<INTMAX_MAX, INTMAX_MAX>::value), INTMAX_MAX );
    BOOST_TEST_EQ( (lcm<INTMAX_MAX, -INTMAX_MAX>::value), INTMAX_MAX );
    BOOST_TEST_EQ( (lcm<-INTMAX_MAX, -INTMAX_MAX>::value), INTMAX_MAX );
    BOOST_TEST_EQ( (lcm<7*11, 11*13>::value), 7*11*13 );

    return boost::report_errors();
}
