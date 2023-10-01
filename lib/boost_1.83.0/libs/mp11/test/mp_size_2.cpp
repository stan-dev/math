// Copyright 2023 Peter Dimov.
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
    using boost::mp11::mp_size;
    using boost::mp11::mp_size_t;

    BOOST_TEST_TRAIT_SAME(mp_size<L1<>>, mp_size_t<0>);
    BOOST_TEST_TRAIT_SAME(mp_size<L1<false>>, mp_size_t<1>);
    BOOST_TEST_TRAIT_SAME(mp_size<L1<false, 0>>, mp_size_t<2>);
    BOOST_TEST_TRAIT_SAME(mp_size<L1<false, 0, std::size_t(0)>>, mp_size_t<3>);

    BOOST_TEST_TRAIT_SAME(mp_size<L2<>>, mp_size_t<0>);
    BOOST_TEST_TRAIT_SAME(mp_size<L2<0>>, mp_size_t<1>);
    BOOST_TEST_TRAIT_SAME(mp_size<L2<0, 1>>, mp_size_t<2>);
    BOOST_TEST_TRAIT_SAME(mp_size<L2<0, 1, 2>>, mp_size_t<3>);

    return boost::report_errors();
}

#endif
