// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/mp11/list.hpp>

#if !defined(BOOST_MP11_HAS_TEMPLATE_AUTO)

#pragma message("Test skipped because BOOST_MP11_HAS_TEMPLATE_AUTO is not defined")
int main() {}

#else

#include <boost/core/lightweight_test_trait.hpp>

template<auto... A> struct L1 {};
template<int... I> struct L2 {};

int main()
{
    using boost::mp11::mp_list_v;
    using boost::mp11::mp_pop_front;
    using boost::mp11::mp_rest;

    BOOST_TEST_TRAIT_SAME(mp_pop_front<mp_list_v<false>>, mp_list_v<>);
    BOOST_TEST_TRAIT_SAME(mp_rest<mp_list_v<false>>, mp_list_v<>);

    BOOST_TEST_TRAIT_SAME(mp_pop_front<L1<-1, true>>, L1<true>);
    BOOST_TEST_TRAIT_SAME(mp_rest<L1<-1, true>>, L1<true>);

    BOOST_TEST_TRAIT_SAME(mp_pop_front<L2<0, 1, 2>>, L2<1, 2>);
    BOOST_TEST_TRAIT_SAME(mp_rest<L2<0, 1, 2>>, L2<1, 2>);

    return boost::report_errors();
}

#endif
