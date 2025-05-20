// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/utility.hpp>

#if !defined(BOOST_MP11_HAS_TEMPLATE_AUTO)

#pragma message("Test skipped because BOOST_MP11_HAS_TEMPLATE_AUTO is not defined")
int main() {}

#else

#include <boost/core/lightweight_test_trait.hpp>

template<auto... A> struct V1 {};
template<int... I> struct V2 {};

int main()
{
    using boost::mp11::mp_pop_back;

    using L1 = V1<false>;
    BOOST_TEST_TRAIT_SAME(mp_pop_back<L1>, V1<>);

    using L2 = V1<false, -1, std::size_t(2)>;
    BOOST_TEST_TRAIT_SAME(mp_pop_back<L2>, V1<false, -1>);

    using L3 = V2<0>;
    BOOST_TEST_TRAIT_SAME(mp_pop_back<L3>, V2<>);

    using L4 = V2<0, 1, 2, 3, 4, 5, 6, 7, 8, 9>;
    BOOST_TEST_TRAIT_SAME(mp_pop_back<L4>, V2<0, 1, 2, 3, 4, 5, 6, 7, 8>);

    using boost::mp11::mp_iota_c;
    using boost::mp11::mp_rename_v;

    int const N = 137;

    using L6 = mp_rename_v<mp_iota_c<N>, V1>;
    using L7 = mp_rename_v<mp_iota_c<N-1>, V1>;
    BOOST_TEST_TRAIT_SAME(mp_pop_back<L6>, L7);

    using boost::mp11::mp_valid;

    using L8 = V2<>;
    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_pop_back, L8>));

    return boost::report_errors();
}

#endif
