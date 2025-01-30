// Copyright 2017, 2023 Peter Dimov.
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
    using boost::mp11::mp_max_element_q;
    using boost::mp11::mp_less;
    using boost::mp11::mp_int;
    using boost::mp11::mp_quote;

    using Q = mp_quote<mp_less>;

    {
        using L1 = V1<1>;
        BOOST_TEST_TRAIT_SAME(mp_max_element_q<L1, Q>, mp_int<1>);

        using L2 = V1<2, 4, 3, 1>;
        BOOST_TEST_TRAIT_SAME(mp_max_element_q<L2, Q>, mp_int<4>);

        using L3 = V1<2, 1, 2, 3>;
        BOOST_TEST_TRAIT_SAME(mp_max_element_q<L3, Q>, mp_int<3>);

        using L4 = V1<-1, 1u, -2, 2u>;
        BOOST_TEST_TRAIT_SAME(mp_max_element_q<L4, Q>, std::integral_constant<unsigned, 2>);
    }

    {
        using L1 = V2<1>;
        BOOST_TEST_TRAIT_SAME(mp_max_element_q<L1, Q>, mp_int<1>);

        using L2 = V2<2, 4, 3, 1>;
        BOOST_TEST_TRAIT_SAME(mp_max_element_q<L2, Q>, mp_int<4>);

        using L3 = V2<2, 1, 2, 3>;
        BOOST_TEST_TRAIT_SAME(mp_max_element_q<L3, Q>, mp_int<3>);

        using L4 = V2<-1, 1, -2, 2>;
        BOOST_TEST_TRAIT_SAME(mp_max_element_q<L4, Q>, mp_int<2>);
    }

    return boost::report_errors();
}

#endif
