// Copyright 2023 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/function.hpp>
#include <boost/core/lightweight_test_trait.hpp>

#if !defined(BOOST_MP11_HAS_TEMPLATE_AUTO)

#pragma message("Test skipped because BOOST_MP11_HAS_TEMPLATE_AUTO is not defined")
int main() {}

#else

template<auto... A> struct V1 {};
template<int... I> struct V2 {};

int main()
{
    using boost::mp11::mp_transform;
    using boost::mp11::mp_identity_t;
    using boost::mp11::mp_plus;

    using L1 = V1<true, -2, 3U, -4L>;
    using L2 = V2<1, -2, 3, -4>;

    BOOST_TEST_TRAIT_SAME(mp_transform<mp_identity_t, L1>, L1);

    BOOST_TEST_TRAIT_SAME(mp_transform<mp_plus, L1>, V1<1, -2, 3U, -4L>);
    BOOST_TEST_TRAIT_SAME(mp_transform<mp_plus, L1, L2>, V1<2, -4, 6U, -8L>);
    BOOST_TEST_TRAIT_SAME(mp_transform<mp_plus, L1, L2, L2>, V1<3, -6, 9U, -12L>);

    BOOST_TEST_TRAIT_SAME(mp_transform<mp_identity_t, L2>, L2);

    BOOST_TEST_TRAIT_SAME(mp_transform<mp_plus, L2>, V2<1, -2, 3, -4>);
    BOOST_TEST_TRAIT_SAME(mp_transform<mp_plus, L2, L1>, V2<2, -4, 6, -8>);
    BOOST_TEST_TRAIT_SAME(mp_transform<mp_plus, L2, L1, L1>, V2<3, -6, 9, -12>);

    //

    return boost::report_errors();
}

#endif
