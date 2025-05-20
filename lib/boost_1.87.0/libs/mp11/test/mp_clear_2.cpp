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
    using boost::mp11::mp_clear;

    BOOST_TEST_TRAIT_SAME(mp_clear<mp_list_v<false, true, 0, -4>>, mp_list_v<>);
    BOOST_TEST_TRAIT_SAME(mp_clear<L1<false, true, 0, -4>>, L1<>);
    BOOST_TEST_TRAIT_SAME(mp_clear<L2<0, 1, 2>>, L2<>);

    return boost::report_errors();
}

#endif
