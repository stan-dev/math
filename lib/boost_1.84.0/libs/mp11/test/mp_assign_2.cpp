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
#include <utility>

template<auto... A> struct L1 {};
template<int... I> struct L2 {};

int main()
{
    using boost::mp11::mp_list;
    using boost::mp11::mp_list_v;
    using boost::mp11::mp_assign;
    using boost::mp11::mp_false;
    using boost::mp11::mp_true;
    using boost::mp11::mp_int;

    //

    BOOST_TEST_TRAIT_SAME(mp_assign<L1<>, mp_list_v<>>, L1<>);
    BOOST_TEST_TRAIT_SAME(mp_assign<L1<>, mp_list_v<true>>, L1<true>);
    BOOST_TEST_TRAIT_SAME(mp_assign<L1<>, mp_list_v<true, -4>>, L1<true, -4>);

    BOOST_TEST_TRAIT_SAME(mp_assign<L1<false, 0>, mp_list_v<>>, L1<>);
    BOOST_TEST_TRAIT_SAME(mp_assign<L1<false, 0>, mp_list_v<false>>, L1<false>);
    BOOST_TEST_TRAIT_SAME(mp_assign<L1<false, 0>, mp_list_v<false, -4>>, L1<false, -4>);

    //

    BOOST_TEST_TRAIT_SAME(mp_assign<L1<>, L2<>>, L1<>);
    BOOST_TEST_TRAIT_SAME(mp_assign<L1<>, L2<0>>, L1<0>);
    BOOST_TEST_TRAIT_SAME(mp_assign<L1<>, L2<0, 1>>, L1<0, 1>);

    BOOST_TEST_TRAIT_SAME(mp_assign<L1<false, 0>, L2<>>, L1<>);
    BOOST_TEST_TRAIT_SAME(mp_assign<L1<false, 0>, L2<0>>, L1<0>);
    BOOST_TEST_TRAIT_SAME(mp_assign<L1<false, 0>, L2<0, 1>>, L1<0, 1>);

    //

    BOOST_TEST_TRAIT_SAME(mp_assign<L2<>, L1<>>, L2<>);
    BOOST_TEST_TRAIT_SAME(mp_assign<L2<>, L1<0>>, L2<0>);
    BOOST_TEST_TRAIT_SAME(mp_assign<L2<>, L1<0, 1>>, L2<0, 1>);

    BOOST_TEST_TRAIT_SAME(mp_assign<L2<1, 2>, L1<>>, L2<>);
    BOOST_TEST_TRAIT_SAME(mp_assign<L2<1, 2>, L1<0>>, L2<0>);
    BOOST_TEST_TRAIT_SAME(mp_assign<L2<1, 2>, L1<0, 1>>, L2<0, 1>);

    //

    BOOST_TEST_TRAIT_SAME(mp_assign<mp_list<>, mp_list_v<>>, mp_list<>);
    BOOST_TEST_TRAIT_SAME(mp_assign<mp_list<>, mp_list_v<true>>, mp_list<mp_true>);
    BOOST_TEST_TRAIT_SAME(mp_assign<mp_list<>, mp_list_v<true, -4>>, mp_list<mp_true, mp_int<-4>>);

    BOOST_TEST_TRAIT_SAME(mp_assign<mp_list<void, float>, mp_list_v<>>, mp_list<>);
    BOOST_TEST_TRAIT_SAME(mp_assign<mp_list<void, float>, mp_list_v<true>>, mp_list<mp_true>);
    BOOST_TEST_TRAIT_SAME(mp_assign<mp_list<void, float>, mp_list_v<true, -4>>, mp_list<mp_true, mp_int<-4>>);

    BOOST_TEST_TRAIT_SAME(mp_assign<std::pair<int, float>, mp_list_v<true, -4>>, std::pair<mp_true, mp_int<-4>>);

    //

    BOOST_TEST_TRAIT_SAME(mp_assign<mp_list_v<>, mp_list<>>, mp_list_v<>);
    BOOST_TEST_TRAIT_SAME(mp_assign<mp_list_v<>, mp_list<mp_true>>, mp_list_v<true>);
    BOOST_TEST_TRAIT_SAME(mp_assign<mp_list_v<>, mp_list<mp_true, mp_int<-4>>>, mp_list_v<true, -4>);

    BOOST_TEST_TRAIT_SAME(mp_assign<mp_list_v<false, 0>, mp_list<>>, mp_list_v<>);
    BOOST_TEST_TRAIT_SAME(mp_assign<mp_list_v<false, 0>, mp_list<mp_true>>, mp_list_v<true>);
    BOOST_TEST_TRAIT_SAME(mp_assign<mp_list_v<false, 0>, mp_list<mp_true, mp_int<-4>>>, mp_list_v<true, -4>);

    BOOST_TEST_TRAIT_SAME(mp_assign<mp_list_v<false, 0>, std::pair<mp_true, mp_int<-4>>>, mp_list_v<true, -4>);

    //

    return boost::report_errors();
}

#endif
