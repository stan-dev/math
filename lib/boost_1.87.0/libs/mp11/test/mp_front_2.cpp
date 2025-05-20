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
    using boost::mp11::mp_list_v;
    using boost::mp11::mp_front;
    using boost::mp11::mp_first;
    using boost::mp11::mp_false;
    using boost::mp11::mp_int;

    BOOST_TEST_TRAIT_SAME(mp_front<mp_list_v<false>>, mp_false);
    BOOST_TEST_TRAIT_SAME(mp_first<mp_list_v<false>>, mp_false);

    BOOST_TEST_TRAIT_SAME(mp_front<L1<-1, true>>, mp_int<-1>);
    BOOST_TEST_TRAIT_SAME(mp_first<L1<-1, true>>, mp_int<-1>);

    BOOST_TEST_TRAIT_SAME(mp_front<L2<0, 1>>, mp_int<0>);
    BOOST_TEST_TRAIT_SAME(mp_first<L2<0, 1>>, mp_int<0>);

    return boost::report_errors();
}

#endif
