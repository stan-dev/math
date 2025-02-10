// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/mp11/list.hpp>

#if !defined(BOOST_MP11_HAS_TEMPLATE_AUTO)

#pragma message("Test skipped because BOOST_MP11_HAS_TEMPLATE_AUTO is not defined")
int main() {}

#else

#include <boost/mp11/integral.hpp>
#include <boost/mp11/utility.hpp>
#include <boost/core/lightweight_test_trait.hpp>

template<auto... A> struct L1 {};
template<int... I> struct L2 {};

template<class N> using double_ = boost::mp11::mp_value< N::value * 2 >;
using Q_double = boost::mp11::mp_quote<double_>;

template<class N> using is_zero = boost::mp11::mp_bool< N::value == 0 >;
using Q_is_zero = boost::mp11::mp_quote<is_zero>;

int main()
{
    using boost::mp11::mp_transform_third;
    using boost::mp11::mp_transform_third_q;

    //

    BOOST_TEST_TRAIT_SAME(mp_transform_third<L1<false, true, 3>, double_>, L1<false, true, 6>);
    BOOST_TEST_TRAIT_SAME(mp_transform_third_q<L1<false, true, 3>, Q_double>, L1<false, true, 6>);

    BOOST_TEST_TRAIT_SAME(mp_transform_third<L1<false, true, 3, 4>, double_>, L1<false, true, 6, 4>);
    BOOST_TEST_TRAIT_SAME(mp_transform_third_q<L1<false, true, 3, 4>, Q_double>, L1<false, true, 6, 4>);

    //

    BOOST_TEST_TRAIT_SAME(mp_transform_third<L1<0, 1, 2>, is_zero>, L1<0, 1, false>);
    BOOST_TEST_TRAIT_SAME(mp_transform_third_q<L1<0, 1, 2>, Q_is_zero>, L1<0, 1, false>);

    BOOST_TEST_TRAIT_SAME(mp_transform_third<L1<0, 1, 2, 3>, is_zero>, L1<0, 1, false, 3>);
    BOOST_TEST_TRAIT_SAME(mp_transform_third_q<L1<0, 1, 2, 3>, Q_is_zero>, L1<0, 1, false, 3>);

    //

    BOOST_TEST_TRAIT_SAME(mp_transform_third<L2<3, 2, 1>, double_>, L2<3, 2, 2>);
    BOOST_TEST_TRAIT_SAME(mp_transform_third_q<L2<3, 2, 1>, Q_double>, L2<3, 2, 2>);

    BOOST_TEST_TRAIT_SAME(mp_transform_third<L2<3, 2, 1, 0>, double_>, L2<3, 2, 2, 0>);
    BOOST_TEST_TRAIT_SAME(mp_transform_third_q<L2<3, 2, 1, 0>, Q_double>, L2<3, 2, 2, 0>);

    //

    return boost::report_errors();
}

#endif
