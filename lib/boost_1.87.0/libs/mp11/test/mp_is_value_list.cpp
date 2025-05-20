// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/mp11/list.hpp>
#include <boost/mp11/integral.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <utility>

template<int... I> struct V2 {};

int main()
{
    using boost::mp11::mp_list;
    using boost::mp11::mp_is_value_list;
    using boost::mp11::mp_true;
    using boost::mp11::mp_false;

    {
        BOOST_TEST_TRAIT_SAME(mp_is_value_list<void>, mp_false);
        BOOST_TEST_TRAIT_SAME(mp_is_value_list<int>, mp_false);
        BOOST_TEST_TRAIT_SAME(mp_is_value_list<char[]>, mp_false);

        BOOST_TEST_TRAIT_SAME(mp_is_value_list<std::pair<void, void>>, mp_false);
    }

    {
        BOOST_TEST_TRAIT_SAME(mp_is_value_list<mp_list<>>, mp_false);
        BOOST_TEST_TRAIT_SAME(mp_is_value_list<mp_list<void>>, mp_false);
        BOOST_TEST_TRAIT_SAME(mp_is_value_list<mp_list<void, void>>, mp_false);
        BOOST_TEST_TRAIT_SAME(mp_is_value_list<mp_list<void, void, void>>, mp_false);
        BOOST_TEST_TRAIT_SAME(mp_is_value_list<mp_list<void, void, void, void>>, mp_false);
    }

#if defined(BOOST_MP11_HAS_TEMPLATE_AUTO)

    using boost::mp11::mp_list_v;

    {
        BOOST_TEST_TRAIT_SAME(mp_is_value_list<mp_list_v<>>, mp_true);
        BOOST_TEST_TRAIT_SAME(mp_is_value_list<mp_list_v<false>>, mp_true);
        BOOST_TEST_TRAIT_SAME(mp_is_value_list<mp_list_v<false, 0>>, mp_true);
        BOOST_TEST_TRAIT_SAME(mp_is_value_list<mp_list_v<false, 0, true>>, mp_true);
        BOOST_TEST_TRAIT_SAME(mp_is_value_list<mp_list_v<false, 0, true, 1>>, mp_true);
    }

    {
        BOOST_TEST_TRAIT_SAME(mp_is_value_list<V2<>>, mp_true);
        BOOST_TEST_TRAIT_SAME(mp_is_value_list<V2<0>>, mp_true);
        BOOST_TEST_TRAIT_SAME(mp_is_value_list<V2<0, 1>>, mp_true);
        BOOST_TEST_TRAIT_SAME(mp_is_value_list<V2<0, 1, 2>>, mp_true);
        BOOST_TEST_TRAIT_SAME(mp_is_value_list<V2<0, 1, 2, 3>>, mp_true);
    }

#endif

    return boost::report_errors();
}
