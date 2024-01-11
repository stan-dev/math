// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/mp11/algorithm.hpp>
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
    using boost::mp11::mp_at;
    using boost::mp11::mp_at_c;
    using boost::mp11::mp_false;
    using boost::mp11::mp_true;
    using boost::mp11::mp_int;
    using boost::mp11::mp_size_t;

    {
        using L1 = V1<false, 0, true, 1, std::size_t(2)>;

        BOOST_TEST_TRAIT_SAME(mp_at_c<L1, 0>, mp_false);
        BOOST_TEST_TRAIT_SAME(mp_at_c<L1, 1>, mp_int<0>);
        BOOST_TEST_TRAIT_SAME(mp_at_c<L1, 2>, mp_true);
        BOOST_TEST_TRAIT_SAME(mp_at_c<L1, 3>, mp_int<1>);
        BOOST_TEST_TRAIT_SAME(mp_at_c<L1, 4>, mp_size_t<2>);

        BOOST_TEST_TRAIT_SAME(mp_at<L1, mp_size_t<0>>, mp_false);
        BOOST_TEST_TRAIT_SAME(mp_at<L1, mp_size_t<1>>, mp_int<0>);
        BOOST_TEST_TRAIT_SAME(mp_at<L1, mp_size_t<2>>, mp_true);
        BOOST_TEST_TRAIT_SAME(mp_at<L1, mp_size_t<3>>, mp_int<1>);
        BOOST_TEST_TRAIT_SAME(mp_at<L1, mp_size_t<4>>, mp_size_t<2>);
    }

    {
        using L1 = V2<1, 2, 3, 4, 5>;

        BOOST_TEST_TRAIT_SAME(mp_at_c<L1, 0>, mp_int<1>);
        BOOST_TEST_TRAIT_SAME(mp_at_c<L1, 1>, mp_int<2>);
        BOOST_TEST_TRAIT_SAME(mp_at_c<L1, 2>, mp_int<3>);
        BOOST_TEST_TRAIT_SAME(mp_at_c<L1, 3>, mp_int<4>);
        BOOST_TEST_TRAIT_SAME(mp_at_c<L1, 4>, mp_int<5>);

        BOOST_TEST_TRAIT_SAME(mp_at<L1, mp_size_t<0>>, mp_int<1>);
        BOOST_TEST_TRAIT_SAME(mp_at<L1, mp_size_t<1>>, mp_int<2>);
        BOOST_TEST_TRAIT_SAME(mp_at<L1, mp_size_t<2>>, mp_int<3>);
        BOOST_TEST_TRAIT_SAME(mp_at<L1, mp_size_t<3>>, mp_int<4>);
        BOOST_TEST_TRAIT_SAME(mp_at<L1, mp_size_t<4>>, mp_int<5>);
    }

    return boost::report_errors();
}

#endif
