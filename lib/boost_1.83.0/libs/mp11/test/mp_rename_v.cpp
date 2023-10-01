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
    using boost::mp11::mp_rename_v;
    using boost::mp11::mp_list_c;
    using boost::mp11::mp_rename;
    using boost::mp11::mp_list;

    BOOST_TEST_TRAIT_SAME(mp_rename_v<L1<>, mp_list_v>, mp_list_v<>);
    BOOST_TEST_TRAIT_SAME(mp_rename_v<L1<false>, mp_list_v>, mp_list_v<false>);
    BOOST_TEST_TRAIT_SAME(mp_rename_v<L1<0>, mp_list_v>, mp_list_v<0>);
    BOOST_TEST_TRAIT_SAME(mp_rename_v<L1<std::size_t(0)>, mp_list_v>, mp_list_v<std::size_t(0)>);

    BOOST_TEST_TRAIT_SAME(mp_rename_v<L1<false, 0, std::size_t(0)>, mp_list_v>, mp_list_v<false, 0, std::size_t(0)>);

    BOOST_TEST_TRAIT_SAME(mp_rename_v<L2<>, mp_list_v>, mp_list_v<>);
    BOOST_TEST_TRAIT_SAME(mp_rename_v<L2<0>, mp_list_v>, mp_list_v<0>);
    BOOST_TEST_TRAIT_SAME(mp_rename_v<L2<0, 1>, mp_list_v>, mp_list_v<0, 1>);
    BOOST_TEST_TRAIT_SAME(mp_rename_v<L2<0, 1, 2>, mp_list_v>, mp_list_v<0, 1, 2>);

    BOOST_TEST_TRAIT_SAME(mp_rename_v<mp_list_v<>, L2>, L2<>);
    BOOST_TEST_TRAIT_SAME(mp_rename_v<mp_list_v<0>, L2>, L2<0>);
    BOOST_TEST_TRAIT_SAME(mp_rename_v<mp_list_v<0, 1>, L2>, L2<0, 1>);
    BOOST_TEST_TRAIT_SAME(mp_rename_v<mp_list_v<0, 1, 2>, L2>, L2<0, 1, 2>);

    BOOST_TEST_TRAIT_SAME(mp_rename_v<mp_list_c<int>, L2>, L2<>);
    BOOST_TEST_TRAIT_SAME(mp_rename_v<mp_list_c<int, 0>, L2>, L2<0>);
    BOOST_TEST_TRAIT_SAME(mp_rename_v<mp_list_c<int, 0, 1>, L2>, L2<0, 1>);
    BOOST_TEST_TRAIT_SAME(mp_rename_v<mp_list_c<int, 0, 1, 2>, L2>, L2<0, 1, 2>);

    using L3 = mp_rename<L1<false, 0, std::size_t(0)>, mp_list>;
    BOOST_TEST_TRAIT_SAME(mp_rename_v<L3, mp_list_v>, mp_list_v<false, 0, std::size_t(0)>);

    return boost::report_errors();
}

#endif
