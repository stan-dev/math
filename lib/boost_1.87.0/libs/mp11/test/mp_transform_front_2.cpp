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
    using boost::mp11::mp_transform_front;
    using boost::mp11::mp_transform_first;
    using boost::mp11::mp_transform_front_q;
    using boost::mp11::mp_transform_first_q;

    //

    BOOST_TEST_TRAIT_SAME(mp_transform_front<L1<3>, double_>, L1<6>);
    BOOST_TEST_TRAIT_SAME(mp_transform_first<L1<3>, double_>, L1<6>);

    BOOST_TEST_TRAIT_SAME(mp_transform_front_q<L1<3>, Q_double>, L1<6>);
    BOOST_TEST_TRAIT_SAME(mp_transform_first_q<L1<3>, Q_double>, L1<6>);

    BOOST_TEST_TRAIT_SAME(mp_transform_front<L1<3, false>, double_>, L1<6, false>);
    BOOST_TEST_TRAIT_SAME(mp_transform_first<L1<3, false>, double_>, L1<6, false>);

    BOOST_TEST_TRAIT_SAME(mp_transform_front_q<L1<3, false>, Q_double>, L1<6, false>);
    BOOST_TEST_TRAIT_SAME(mp_transform_first_q<L1<3, false>, Q_double>, L1<6, false>);

    BOOST_TEST_TRAIT_SAME(mp_transform_front<L1<3, false, true>, double_>, L1<6, false, true>);
    BOOST_TEST_TRAIT_SAME(mp_transform_first<L1<3, false, true>, double_>, L1<6, false, true>);

    BOOST_TEST_TRAIT_SAME(mp_transform_front_q<L1<3, false, true>, Q_double>, L1<6, false, true>);
    BOOST_TEST_TRAIT_SAME(mp_transform_first_q<L1<3, false, true>, Q_double>, L1<6, false, true>);

    //

    BOOST_TEST_TRAIT_SAME(mp_transform_front<L1<0>, is_zero>, L1<true>);
    BOOST_TEST_TRAIT_SAME(mp_transform_first<L1<0>, is_zero>, L1<true>);

    BOOST_TEST_TRAIT_SAME(mp_transform_front_q<L1<0>, Q_is_zero>, L1<true>);
    BOOST_TEST_TRAIT_SAME(mp_transform_first_q<L1<0>, Q_is_zero>, L1<true>);

    BOOST_TEST_TRAIT_SAME(mp_transform_front<L1<0, 1>, is_zero>, L1<true, 1>);
    BOOST_TEST_TRAIT_SAME(mp_transform_first<L1<0, 1>, is_zero>, L1<true, 1>);

    BOOST_TEST_TRAIT_SAME(mp_transform_front_q<L1<0, 1>, Q_is_zero>, L1<true, 1>);
    BOOST_TEST_TRAIT_SAME(mp_transform_first_q<L1<0, 1>, Q_is_zero>, L1<true, 1>);

    BOOST_TEST_TRAIT_SAME(mp_transform_front<L1<0, 1, 2>, is_zero>, L1<true, 1, 2>);
    BOOST_TEST_TRAIT_SAME(mp_transform_first<L1<0, 1, 2>, is_zero>, L1<true, 1, 2>);

    BOOST_TEST_TRAIT_SAME(mp_transform_front_q<L1<0, 1, 2>, Q_is_zero>, L1<true, 1, 2>);
    BOOST_TEST_TRAIT_SAME(mp_transform_first_q<L1<0, 1, 2>, Q_is_zero>, L1<true, 1, 2>);

    //

    BOOST_TEST_TRAIT_SAME(mp_transform_front<L2<3>, double_>, L2<6>);
    BOOST_TEST_TRAIT_SAME(mp_transform_first<L2<3>, double_>, L2<6>);

    BOOST_TEST_TRAIT_SAME(mp_transform_front_q<L2<3>, Q_double>, L2<6>);
    BOOST_TEST_TRAIT_SAME(mp_transform_first_q<L2<3>, Q_double>, L2<6>);

    BOOST_TEST_TRAIT_SAME(mp_transform_front<L2<3, 2>, double_>, L2<6, 2>);
    BOOST_TEST_TRAIT_SAME(mp_transform_first<L2<3, 2>, double_>, L2<6, 2>);

    BOOST_TEST_TRAIT_SAME(mp_transform_front_q<L2<3, 2>, Q_double>, L2<6, 2>);
    BOOST_TEST_TRAIT_SAME(mp_transform_first_q<L2<3, 2>, Q_double>, L2<6, 2>);

    BOOST_TEST_TRAIT_SAME(mp_transform_front<L2<3, 2, 1>, double_>, L2<6, 2, 1>);
    BOOST_TEST_TRAIT_SAME(mp_transform_first<L2<3, 2, 1>, double_>, L2<6, 2, 1>);

    BOOST_TEST_TRAIT_SAME(mp_transform_front_q<L2<3, 2, 1>, Q_double>, L2<6, 2, 1>);
    BOOST_TEST_TRAIT_SAME(mp_transform_first_q<L2<3, 2, 1>, Q_double>, L2<6, 2, 1>);

    //

    return boost::report_errors();
}

#endif
