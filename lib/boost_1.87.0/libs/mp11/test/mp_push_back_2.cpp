// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/mp11/list.hpp>

#if !defined(BOOST_MP11_HAS_TEMPLATE_AUTO)

#pragma message("Test skipped because BOOST_MP11_HAS_TEMPLATE_AUTO is not defined")
int main() {}

#else

#include <boost/mp11/integral.hpp>
#include <boost/core/lightweight_test_trait.hpp>

template<auto... A> struct L1 {};
template<int... I> struct L2 {};

int main()
{
    using boost::mp11::mp_push_back;
    using boost::mp11::mp_true;
    using boost::mp11::mp_false;
    using boost::mp11::mp_int;

    BOOST_TEST_TRAIT_SAME(mp_push_back<L1<>>, L1<>);
    BOOST_TEST_TRAIT_SAME(mp_push_back<L1<>, mp_false>, L1<false>);
    BOOST_TEST_TRAIT_SAME(mp_push_back<L1<>, mp_false, mp_true>, L1<false, true>);

    BOOST_TEST_TRAIT_SAME(mp_push_back<L1<-1, -2>>, L1<-1, -2>);
    BOOST_TEST_TRAIT_SAME(mp_push_back<L1<-1, -2>, mp_false>, L1<-1, -2, false>);
    BOOST_TEST_TRAIT_SAME(mp_push_back<L1<-1, -2>, mp_false, mp_true>, L1<-1, -2, false, true>);

    BOOST_TEST_TRAIT_SAME(mp_push_back<L2<>>, L2<>);
    BOOST_TEST_TRAIT_SAME(mp_push_back<L2<>, mp_int<0>>, L2<0>);
    BOOST_TEST_TRAIT_SAME(mp_push_back<L2<>, mp_int<0>, mp_int<1>>, L2<0, 1>);

    BOOST_TEST_TRAIT_SAME(mp_push_back<L2<-1, -2>>, L2<-1, -2>);
    BOOST_TEST_TRAIT_SAME(mp_push_back<L2<-1, -2>, mp_int<0>>, L2<-1, -2, 0>);
    BOOST_TEST_TRAIT_SAME(mp_push_back<L2<-1, -2>, mp_int<0>, mp_int<1>>, L2<-1, -2, 0, 1>);

    return boost::report_errors();
}

#endif
