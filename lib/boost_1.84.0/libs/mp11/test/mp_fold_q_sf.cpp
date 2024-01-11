// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/integral.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <type_traits>
#include <tuple>

using boost::mp11::mp_size_t;

struct Q
{
    template<class N, class T> using fn = mp_size_t<N::value + sizeof(T)>;
};

int main()
{
    using boost::mp11::mp_valid;
    using boost::mp11::mp_fold_q;
    using boost::mp11::mp_list;

    BOOST_TEST_TRAIT_TRUE((mp_valid<mp_fold_q, mp_list<>, mp_size_t<0>, Q>));
    BOOST_TEST_TRAIT_TRUE((mp_valid<mp_fold_q, mp_list<int>, mp_size_t<0>, Q>));
    BOOST_TEST_TRAIT_TRUE((mp_valid<mp_fold_q, mp_list<int, int>, mp_size_t<0>, Q>));
    BOOST_TEST_TRAIT_TRUE((mp_valid<mp_fold_q, mp_list<int, int, int>, mp_size_t<0>, Q>));
    BOOST_TEST_TRAIT_TRUE((mp_valid<mp_fold_q, mp_list<int, int, int, int>, mp_size_t<0>, Q>));

    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_fold_q, mp_list<void>, mp_size_t<0>, Q>));
    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_fold_q, mp_list<int, void>, mp_size_t<0>, Q>));
    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_fold_q, mp_list<int, int, void>, mp_size_t<0>, Q>));
    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_fold_q, mp_list<int, int, int, void>, mp_size_t<0>, Q>));

#if !BOOST_MP11_WORKAROUND( BOOST_MP11_MSVC, < 1910 )

    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_fold_q, mp_list<void, int>, mp_size_t<0>, Q>));
    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_fold_q, mp_list<void, int, int>, mp_size_t<0>, Q>));
    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_fold_q, mp_list<void, int, int, int>, mp_size_t<0>, Q>));

    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_fold_q, mp_list<int, void, int>, mp_size_t<0>, Q>));
    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_fold_q, mp_list<int, void, int, int>, mp_size_t<0>, Q>));

    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_fold_q, mp_list<int, int, void, int>, mp_size_t<0>, Q>));

#endif

    return boost::report_errors();
}
