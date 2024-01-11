// Copyright 2023 Peter Dimov
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
    using boost::mp11::mp_insert;
    using boost::mp11::mp_insert_c;
    using boost::mp11::mp_size_t;
    using boost::mp11::mp_false;
    using boost::mp11::mp_true;
    using boost::mp11::mp_int;

    {
        using L1 = V1<>;

        BOOST_TEST_TRAIT_SAME(mp_insert_c<L1, 0>, L1);
        BOOST_TEST_TRAIT_SAME(mp_insert<L1, mp_size_t<0>>, L1);

        BOOST_TEST_TRAIT_SAME(mp_insert_c<L1, 0, mp_false>, V1<false>);
        BOOST_TEST_TRAIT_SAME(mp_insert<L1, mp_size_t<0>, mp_false>, V1<false>);

        BOOST_TEST_TRAIT_SAME(mp_insert_c<L1, 0, mp_false, mp_true>, V1<false, true>);
        BOOST_TEST_TRAIT_SAME(mp_insert<L1, mp_size_t<0>, mp_false, mp_true>, V1<false, true>);

        using L2 = V1<1, 2, 3, 4, 5>;

        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 0>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 1>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 2>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 3>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 4>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 5>, L2);

        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 0, mp_false, mp_true>, V1<false, true, 1, 2, 3, 4, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 1, mp_false, mp_true>, V1<1, false, true, 2, 3, 4, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 2, mp_false, mp_true>, V1<1, 2, false, true, 3, 4, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 3, mp_false, mp_true>, V1<1, 2, 3, false, true, 4, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 4, mp_false, mp_true>, V1<1, 2, 3, 4, false, true, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 5, mp_false, mp_true>, V1<1, 2, 3, 4, 5, false, true>);

        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<0>>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<1>>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<2>>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<3>>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<4>>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<5>>, L2);

        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<0>, mp_false, mp_true>, V1<false, true, 1, 2, 3, 4, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<1>, mp_false, mp_true>, V1<1, false, true, 2, 3, 4, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<2>, mp_false, mp_true>, V1<1, 2, false, true, 3, 4, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<3>, mp_false, mp_true>, V1<1, 2, 3, false, true, 4, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<4>, mp_false, mp_true>, V1<1, 2, 3, 4, false, true, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<5>, mp_false, mp_true>, V1<1, 2, 3, 4, 5, false, true>);
    }

    {
        using L1 = V2<>;

        BOOST_TEST_TRAIT_SAME(mp_insert_c<L1, 0>, L1);
        BOOST_TEST_TRAIT_SAME(mp_insert<L1, mp_size_t<0>>, L1);

        BOOST_TEST_TRAIT_SAME(mp_insert_c<L1, 0, mp_int<-1>>, V2<-1>);
        BOOST_TEST_TRAIT_SAME(mp_insert<L1, mp_size_t<0>, mp_int<-1>>, V2<-1>);

        BOOST_TEST_TRAIT_SAME(mp_insert_c<L1, 0, mp_int<-1>, mp_int<-2>>, V2<-1, -2>);
        BOOST_TEST_TRAIT_SAME(mp_insert<L1, mp_size_t<0>, mp_int<-1>, mp_int<-2>>, V2<-1, -2>);

        using L2 = V2<1, 2, 3, 4, 5>;

        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 0>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 1>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 2>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 3>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 4>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 5>, L2);

        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 0, mp_int<-1>, mp_int<-2>>, V2<-1, -2, 1, 2, 3, 4, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 1, mp_int<-1>, mp_int<-2>>, V2<1, -1, -2, 2, 3, 4, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 2, mp_int<-1>, mp_int<-2>>, V2<1, 2, -1, -2, 3, 4, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 3, mp_int<-1>, mp_int<-2>>, V2<1, 2, 3, -1, -2, 4, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 4, mp_int<-1>, mp_int<-2>>, V2<1, 2, 3, 4, -1, -2, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert_c<L2, 5, mp_int<-1>, mp_int<-2>>, V2<1, 2, 3, 4, 5, -1, -2>);

        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<0>>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<1>>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<2>>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<3>>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<4>>, L2);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<5>>, L2);

        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<0>, mp_int<-1>, mp_int<-2>>, V2<-1, -2, 1, 2, 3, 4, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<1>, mp_int<-1>, mp_int<-2>>, V2<1, -1, -2, 2, 3, 4, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<2>, mp_int<-1>, mp_int<-2>>, V2<1, 2, -1, -2, 3, 4, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<3>, mp_int<-1>, mp_int<-2>>, V2<1, 2, 3, -1, -2, 4, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<4>, mp_int<-1>, mp_int<-2>>, V2<1, 2, 3, 4, -1, -2, 5>);
        BOOST_TEST_TRAIT_SAME(mp_insert<L2, mp_size_t<5>, mp_int<-1>, mp_int<-2>>, V2<1, 2, 3, 4, 5, -1, -2>);
    }

    return boost::report_errors();
}

#endif
