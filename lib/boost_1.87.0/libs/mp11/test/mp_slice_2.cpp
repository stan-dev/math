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
    using boost::mp11::mp_slice;
    using boost::mp11::mp_slice_c;
    using boost::mp11::mp_size_t;

    {
        using L1 = V1<>;

        BOOST_TEST_TRAIT_SAME(mp_slice_c<L1, 0, 0>, L1);
        BOOST_TEST_TRAIT_SAME(mp_slice<L1, mp_size_t<0>, mp_size_t<0>>, L1);

        using L2 = V1<false, 0, true, 1, std::size_t(2)>;

        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 0, 0>, V1<>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 0, 1>, V1<false>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 0, 2>, V1<false, 0>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 0, 3>, V1<false, 0, true>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 0, 4>, V1<false, 0, true, 1>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 0, 5>, V1<false, 0, true, 1, std::size_t(2)>);

        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 1, 1>, V1<>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 1, 2>, V1<0>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 1, 3>, V1<0, true>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 1, 4>, V1<0, true, 1>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 1, 5>, V1<0, true, 1, std::size_t(2)>);

        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 2, 2>, V1<>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 2, 3>, V1<true>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 2, 4>, V1<true, 1>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 2, 5>, V1<true, 1, std::size_t(2)>);

        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 3, 3>, V1<>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 3, 4>, V1<1>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 3, 5>, V1<1, std::size_t(2)>);

        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 4, 4>, V1<>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 4, 5>, V1<std::size_t(2)>);

        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 5, 5>, V1<>);

        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<0>, mp_size_t<0>>, V1<>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<0>, mp_size_t<1>>, V1<false>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<0>, mp_size_t<2>>, V1<false, 0>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<0>, mp_size_t<3>>, V1<false, 0, true>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<0>, mp_size_t<4>>, V1<false, 0, true, 1>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<0>, mp_size_t<5>>, V1<false, 0, true, 1, std::size_t(2)>);

        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<1>, mp_size_t<1>>, V1<>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<1>, mp_size_t<2>>, V1<0>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<1>, mp_size_t<3>>, V1<0, true>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<1>, mp_size_t<4>>, V1<0, true, 1>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<1>, mp_size_t<5>>, V1<0, true, 1, std::size_t(2)>);

        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<2>, mp_size_t<2>>, V1<>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<2>, mp_size_t<3>>, V1<true>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<2>, mp_size_t<4>>, V1<true, 1>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<2>, mp_size_t<5>>, V1<true, 1, std::size_t(2)>);

        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<3>, mp_size_t<3>>, V1<>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<3>, mp_size_t<4>>, V1<1>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<3>, mp_size_t<5>>, V1<1, std::size_t(2)>);

        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<4>, mp_size_t<4>>, V1<>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<4>, mp_size_t<5>>, V1<std::size_t(2)>);

        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<5>, mp_size_t<5>>, V1<>);
    }

    {
        using L1 = V2<>;

        BOOST_TEST_TRAIT_SAME(mp_slice_c<L1, 0, 0>, L1);
        BOOST_TEST_TRAIT_SAME(mp_slice<L1, mp_size_t<0>, mp_size_t<0>>, L1);

        using L2 = V2<1, 2, 3, 4, 5>;

        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 0, 0>, V2<>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 0, 1>, V2<1>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 0, 2>, V2<1, 2>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 0, 3>, V2<1, 2, 3>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 0, 4>, V2<1, 2, 3, 4>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 0, 5>, V2<1, 2, 3, 4, 5>);

        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 1, 1>, V2<>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 1, 2>, V2<2>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 1, 3>, V2<2, 3>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 1, 4>, V2<2, 3, 4>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 1, 5>, V2<2, 3, 4, 5>);

        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 2, 2>, V2<>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 2, 3>, V2<3>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 2, 4>, V2<3, 4>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 2, 5>, V2<3, 4, 5>);

        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 3, 3>, V2<>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 3, 4>, V2<4>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 3, 5>, V2<4, 5>);

        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 4, 4>, V2<>);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 4, 5>, V2<5>);

        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 5, 5>, V2<>);

        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<0>, mp_size_t<0>>, V2<>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<0>, mp_size_t<1>>, V2<1>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<0>, mp_size_t<2>>, V2<1, 2>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<0>, mp_size_t<3>>, V2<1, 2, 3>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<0>, mp_size_t<4>>, V2<1, 2, 3, 4>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<0>, mp_size_t<5>>, V2<1, 2, 3, 4, 5>);

        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<1>, mp_size_t<1>>, V2<>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<1>, mp_size_t<2>>, V2<2>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<1>, mp_size_t<3>>, V2<2, 3>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<1>, mp_size_t<4>>, V2<2, 3, 4>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<1>, mp_size_t<5>>, V2<2, 3, 4, 5>);

        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<2>, mp_size_t<2>>, V2<>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<2>, mp_size_t<3>>, V2<3>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<2>, mp_size_t<4>>, V2<3, 4>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<2>, mp_size_t<5>>, V2<3, 4, 5>);

        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<3>, mp_size_t<3>>, V2<>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<3>, mp_size_t<4>>, V2<4>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<3>, mp_size_t<5>>, V2<4, 5>);

        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<4>, mp_size_t<4>>, V2<>);
        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<4>, mp_size_t<5>>, V2<5>);

        BOOST_TEST_TRAIT_SAME(mp_slice<L2, mp_size_t<5>, mp_size_t<5>>, V2<>);
    }

    using boost::mp11::mp_iota_c;
    using boost::mp11::mp_rename_v;

    {
        using L1 = mp_rename_v<mp_iota_c<71>, V1>;
        using L2 = mp_rename_v<mp_iota_c<72>, V1>;

        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 0, 72>, L2);
        BOOST_TEST_TRAIT_SAME(mp_slice_c<L2, 0, 71>, L1);
    }

    return boost::report_errors();
}

#endif
