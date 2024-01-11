
// Copyright 2017 Peter Dimov.
//
// Distributed under the Boost Software License, Version 1.0.
//
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt


#include <boost/mp11/list.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <type_traits>
#include <tuple>
#include <utility>

template<int... I> struct V2 {};

int main()
{
    using boost::mp11::mp_list;
    using boost::mp11::mp_is_list;
    using boost::mp11::mp_true;
    using boost::mp11::mp_false;

    {
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_is_list<void>, mp_false>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_is_list<int>, mp_false>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_is_list<char[]>, mp_false>));
    }

    {
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_is_list<mp_list<>>, mp_true>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_is_list<mp_list<void>>, mp_true>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_is_list<mp_list<void, void>>, mp_true>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_is_list<mp_list<void, void, void>>, mp_true>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_is_list<mp_list<void, void, void, void>>, mp_true>));
    }

    {
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_is_list<std::tuple<>>, mp_true>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_is_list<std::tuple<void>>, mp_true>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_is_list<std::tuple<void, void>>, mp_true>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_is_list<std::tuple<void, void, void>>, mp_true>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_is_list<std::tuple<void, void, void, void>>, mp_true>));
    }

    {
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_is_list<std::pair<void, void>>, mp_true>));
    }

#if defined(BOOST_MP11_HAS_TEMPLATE_AUTO)

    using boost::mp11::mp_list_v;

    {
        BOOST_TEST_TRAIT_SAME(mp_is_list<mp_list_v<>>, mp_false);
        BOOST_TEST_TRAIT_SAME(mp_is_list<mp_list_v<false>>, mp_false);
        BOOST_TEST_TRAIT_SAME(mp_is_list<mp_list_v<false, 0>>, mp_false);
        BOOST_TEST_TRAIT_SAME(mp_is_list<mp_list_v<false, 0, true>>, mp_false);
        BOOST_TEST_TRAIT_SAME(mp_is_list<mp_list_v<false, 0, true, 1>>, mp_false);
    }

    {
        BOOST_TEST_TRAIT_SAME(mp_is_list<V2<>>, mp_false);
        BOOST_TEST_TRAIT_SAME(mp_is_list<V2<0>>, mp_false);
        BOOST_TEST_TRAIT_SAME(mp_is_list<V2<0, 1>>, mp_false);
        BOOST_TEST_TRAIT_SAME(mp_is_list<V2<0, 1, 2>>, mp_false);
        BOOST_TEST_TRAIT_SAME(mp_is_list<V2<0, 1, 2, 3>>, mp_false);
    }

#endif

    return boost::report_errors();
}
