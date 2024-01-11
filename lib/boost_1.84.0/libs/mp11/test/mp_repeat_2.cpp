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
    using boost::mp11::mp_repeat;
    using boost::mp11::mp_repeat_c;
    using boost::mp11::mp_false;
    using boost::mp11::mp_true;
    using boost::mp11::mp_int;
    using boost::mp11::mp_size_t;

    using L1 = V1<>;

    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L1, 0>, V1<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L1, 1>, V1<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L1, 2>, V1<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L1, 3>, V1<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L1, 31>, V1<>);

    BOOST_TEST_TRAIT_SAME(mp_repeat<L1, mp_false>, V1<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L1, mp_true>, V1<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L1, mp_int<2>>, V1<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L1, mp_int<3>>, V1<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L1, mp_size_t<31>>, V1<>);

    using L2 = V1<true>;

    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L2, 0>, V1<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L2, 1>, V1<true>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L2, 2>, V1<true, true>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L2, 3>, V1<true, true, true>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L2, 4>, V1<true, true, true, true>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L2, 5>, V1<true, true, true, true, true>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L2, 6>, V1<true, true, true, true, true, true>);

    BOOST_TEST_TRAIT_SAME(mp_repeat<L2, mp_false>, V1<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L2, mp_true>, V1<true>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L2, mp_int<2>>, V1<true, true>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L2, mp_int<3>>, V1<true, true, true>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L2, mp_int<4>>, V1<true, true, true, true>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L2, mp_size_t<5>>, V1<true, true, true, true, true>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L2, mp_size_t<6>>, V1<true, true, true, true, true, true>);

    using L3 = V1<true, 2>;

    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L3, 0>, V1<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L3, 1>, V1<true, 2>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L3, 2>, V1<true, 2, true, 2>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L3, 3>, V1<true, 2, true, 2, true, 2>);

    BOOST_TEST_TRAIT_SAME(mp_repeat<L3, mp_false>, V1<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L3, mp_true>, V1<true, 2>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L3, mp_int<2>>, V1<true, 2, true, 2>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L3, mp_size_t<3>>, V1<true, 2, true, 2, true, 2>);

    //

    using L4 = V2<>;

    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L4, 0>, V2<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L4, 1>, V2<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L4, 2>, V2<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L4, 3>, V2<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L4, 31>, V2<>);

    BOOST_TEST_TRAIT_SAME(mp_repeat<L4, mp_false>, V2<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L4, mp_true>, V2<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L4, mp_int<2>>, V2<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L4, mp_int<3>>, V2<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L4, mp_size_t<31>>, V2<>);

    using L5 = V2<1>;

    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L5, 0>, V2<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L5, 1>, V2<1>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L5, 2>, V2<1, 1>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L5, 3>, V2<1, 1, 1>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L5, 4>, V2<1, 1, 1, 1>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L5, 5>, V2<1, 1, 1, 1, 1>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L5, 6>, V2<1, 1, 1, 1, 1, 1>);

    BOOST_TEST_TRAIT_SAME(mp_repeat<L5, mp_false>, V2<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L5, mp_true>, V2<1>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L5, mp_int<2>>, V2<1, 1>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L5, mp_int<3>>, V2<1, 1, 1>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L5, mp_int<4>>, V2<1, 1, 1, 1>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L5, mp_size_t<5>>, V2<1, 1, 1, 1, 1>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L5, mp_size_t<6>>, V2<1, 1, 1, 1, 1, 1>);

    using L6 = V2<1, 2>;

    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L6, 0>, V2<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L6, 1>, V2<1, 2>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L6, 2>, V2<1, 2, 1, 2>);
    BOOST_TEST_TRAIT_SAME(mp_repeat_c<L6, 3>, V2<1, 2, 1, 2, 1, 2>);

    BOOST_TEST_TRAIT_SAME(mp_repeat<L6, mp_false>, V2<>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L6, mp_true>, V2<1, 2>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L6, mp_int<2>>, V2<1, 2, 1, 2>);
    BOOST_TEST_TRAIT_SAME(mp_repeat<L6, mp_size_t<3>>, V2<1, 2, 1, 2, 1, 2>);

    //

    return boost::report_errors();
}

#endif
