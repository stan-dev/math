// Copyright 2015-2023 Peter Dimov.
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
    using boost::mp11::mp_list;
    using boost::mp11::mp_rename;
    using boost::mp11::mp_true;
    using boost::mp11::mp_false;
    using boost::mp11::mp_int;
    using boost::mp11::mp_size_t;
    using boost::mp11::mp_list_c;

    BOOST_TEST_TRAIT_SAME(mp_rename<L1<>, mp_list>, mp_list<>);

    BOOST_TEST_TRAIT_SAME(mp_rename<L1<false>, mp_list>, mp_list<mp_false>);
    BOOST_TEST_TRAIT_SAME(mp_rename<L1<true>, mp_list>, mp_list<mp_true>);

    BOOST_TEST_TRAIT_SAME(mp_rename<L1<0>, mp_list>, mp_list<mp_int<0>>);
    BOOST_TEST_TRAIT_SAME(mp_rename<L1<-1>, mp_list>, mp_list<mp_int<-1>>);

    BOOST_TEST_TRAIT_SAME(mp_rename<L1<std::size_t(0)>, mp_list>, mp_list<mp_size_t<0>>);
    BOOST_TEST_TRAIT_SAME(mp_rename<L1<std::size_t(1)>, mp_list>, mp_list<mp_size_t<1>>);

    BOOST_TEST_TRAIT_SAME(mp_rename<L1<false, 0, std::size_t(0)>, mp_list>, mp_list<mp_false, mp_int<0>, mp_size_t<0>>);

    BOOST_TEST_TRAIT_SAME(mp_rename<L2<>, mp_list>, mp_list_c<int>);
    BOOST_TEST_TRAIT_SAME(mp_rename<L2<0>, mp_list>, mp_list_c<int, 0>);
    BOOST_TEST_TRAIT_SAME(mp_rename<L2<0, 1>, mp_list>, mp_list_c<int, 0, 1>);
    BOOST_TEST_TRAIT_SAME(mp_rename<L2<0, 1, 2>, mp_list>, mp_list_c<int, 0, 1, 2>);

    return boost::report_errors();
}

#endif
