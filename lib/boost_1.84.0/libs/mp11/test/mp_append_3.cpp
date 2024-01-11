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
    using boost::mp11::mp_append;

    //

    {
        using L1 = mp_list_v<false, 1>;
        using L2 = mp_list_v<false, 2>;
        using L3 = mp_list_v<false, 3>;
        using L4 = mp_list_v<false, 4>;

        BOOST_TEST_TRAIT_SAME(mp_append<L1>, mp_list_v<false, 1>);
        BOOST_TEST_TRAIT_SAME(mp_append<L1, L2>, mp_list_v<false, 1, false, 2>);
        BOOST_TEST_TRAIT_SAME(mp_append<L1, L2, L3>, mp_list_v<false, 1, false, 2, false, 3>);
        BOOST_TEST_TRAIT_SAME(mp_append<L1, L2, L3, L4>, mp_list_v<false, 1, false, 2, false, 3, false, 4>);
    }

    //

    BOOST_TEST_TRAIT_SAME(mp_append<L1<>>, L1<>);
    BOOST_TEST_TRAIT_SAME(mp_append<L1<>, L2<1>>, L1<1>);
    BOOST_TEST_TRAIT_SAME(mp_append<L1<>, L2<1>, L1<true>>, L1<1, true>);
    BOOST_TEST_TRAIT_SAME(mp_append<L1<>, L2<1>, L1<true>, L2<2>>, L1<1, true, 2>);
    BOOST_TEST_TRAIT_SAME(mp_append<L1<>, L2<1>, L1<true>, L2<2>, L1<false>>, L1<1, true, 2, false>);
    BOOST_TEST_TRAIT_SAME(mp_append<L1<>, L2<1>, L1<true>, L2<2>, L1<false>, L2<3>>, L1<1, true, 2, false, 3>);

    //

    BOOST_TEST_TRAIT_SAME(mp_append<L2<1>>, L2<1>);
    BOOST_TEST_TRAIT_SAME(mp_append<L2<1>, L2<2>>, L2<1, 2>);
    BOOST_TEST_TRAIT_SAME(mp_append<L2<1>, L2<2>, L2<3>>, L2<1, 2, 3>);
    BOOST_TEST_TRAIT_SAME(mp_append<L2<1>, L2<2>, L2<3>, L2<4>>, L2<1, 2, 3, 4>);
    BOOST_TEST_TRAIT_SAME(mp_append<L2<1>, L2<2>, L2<3>, L2<4>, L2<5>>, L2<1, 2, 3, 4, 5>);

    //

    return boost::report_errors();
}

#endif
