// Copyright 2023 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/integral.hpp>

#if !defined(BOOST_MP11_HAS_TEMPLATE_AUTO)

#pragma message("Test skipped because BOOST_MP11_HAS_TEMPLATE_AUTO is not defined")
int main() {}

#else

#include <boost/core/lightweight_test_trait.hpp>

template<auto... A> struct V1 {};
template<int... I> struct V2 {};

int main()
{
    using boost::mp11::mp_list;
    using boost::mp11::mp_take;
    using boost::mp11::mp_take_c;
    using boost::mp11::mp_size_t;

    {
        using L1 = V1<>;

        BOOST_TEST_TRAIT_SAME(mp_take_c<L1, 0>, L1);
        BOOST_TEST_TRAIT_SAME(mp_take<L1, mp_size_t<0>>, L1);

        using L2 = V1<false, 0, true, 1, std::size_t(2)>;

        BOOST_TEST_TRAIT_SAME(mp_take_c<L2, 0>, V1<>);
        BOOST_TEST_TRAIT_SAME(mp_take_c<L2, 1>, V1<false>);
        BOOST_TEST_TRAIT_SAME(mp_take_c<L2, 2>, V1<false, 0>);
        BOOST_TEST_TRAIT_SAME(mp_take_c<L2, 3>, V1<false, 0, true>);
        BOOST_TEST_TRAIT_SAME(mp_take_c<L2, 4>, V1<false, 0, true, 1>);
        BOOST_TEST_TRAIT_SAME(mp_take_c<L2, 5>, V1<false, 0, true, 1, std::size_t(2)>);

        BOOST_TEST_TRAIT_SAME(mp_take<L2, mp_size_t<0>>, V1<>);
        BOOST_TEST_TRAIT_SAME(mp_take<L2, mp_size_t<1>>, V1<false>);
        BOOST_TEST_TRAIT_SAME(mp_take<L2, mp_size_t<2>>, V1<false, 0>);
        BOOST_TEST_TRAIT_SAME(mp_take<L2, mp_size_t<3>>, V1<false, 0, true>);
        BOOST_TEST_TRAIT_SAME(mp_take<L2, mp_size_t<4>>, V1<false, 0, true, 1>);
        BOOST_TEST_TRAIT_SAME(mp_take<L2, mp_size_t<5>>, V1<false, 0, true, 1, std::size_t(2)>);
    }

    {
        using L1 = V2<>;

        BOOST_TEST_TRAIT_SAME(mp_take_c<L1, 0>, L1);
        BOOST_TEST_TRAIT_SAME(mp_take<L1, mp_size_t<0>>, L1);

        using L2 = V2<1, 2, 3, 4, 5>;

        BOOST_TEST_TRAIT_SAME(mp_take_c<L2, 0>, V2<>);
        BOOST_TEST_TRAIT_SAME(mp_take_c<L2, 1>, V2<1>);
        BOOST_TEST_TRAIT_SAME(mp_take_c<L2, 2>, V2<1, 2>);
        BOOST_TEST_TRAIT_SAME(mp_take_c<L2, 3>, V2<1, 2, 3>);
        BOOST_TEST_TRAIT_SAME(mp_take_c<L2, 4>, V2<1, 2, 3, 4>);
        BOOST_TEST_TRAIT_SAME(mp_take_c<L2, 5>, V2<1, 2, 3, 4, 5>);

        BOOST_TEST_TRAIT_SAME(mp_take<L2, mp_size_t<0>>, V2<>);
        BOOST_TEST_TRAIT_SAME(mp_take<L2, mp_size_t<1>>, V2<1>);
        BOOST_TEST_TRAIT_SAME(mp_take<L2, mp_size_t<2>>, V2<1, 2>);
        BOOST_TEST_TRAIT_SAME(mp_take<L2, mp_size_t<3>>, V2<1, 2, 3>);
        BOOST_TEST_TRAIT_SAME(mp_take<L2, mp_size_t<4>>, V2<1, 2, 3, 4>);
        BOOST_TEST_TRAIT_SAME(mp_take<L2, mp_size_t<5>>, V2<1, 2, 3, 4, 5>);
    }

    using boost::mp11::mp_iota_c;
    using boost::mp11::mp_rename_v;

    {
        using L1 = mp_rename_v<mp_iota_c<71>, V1>;
        using L2 = mp_rename_v<mp_iota_c<72>, V1>;

        BOOST_TEST_TRAIT_SAME(mp_take_c<L2, 72>, L2);
        BOOST_TEST_TRAIT_SAME(mp_take_c<L2, 71>, L1);
    }

    return boost::report_errors();
}

#endif
