// Copyright 2023 Braden Ganetsky (braden.ganetsky@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0.
//
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/integral.hpp>
#include <boost/mp11/function.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <tuple>
#include <utility>

struct X1 {};
struct X2 {};
struct X3 {};
struct X4 {};
struct X5 {};

using boost::mp11::mp_plus;
using boost::mp11::mp_int;

#if !BOOST_MP11_WORKAROUND( BOOST_MP11_MSVC, <= 1800 )

struct average
{
    template<class... C> using fn = mp_int<mp_plus<C...>::value / sizeof...(C)>;
};

#else

template<class... C> struct average_impl: mp_int<mp_plus<C...>::value / sizeof...(C)> {};

struct average
{
    template<class... C> using fn = typename average_impl<C...>::type;
};

#endif

int main()
{
    using boost::mp11::mp_list;
    using boost::mp11::mp_list_c;
    using boost::mp11::mp_quote;
    using boost::mp11::mp_size_t;
    using boost::mp11::mp_sliding_fold_q;

    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<>, mp_size_t<1>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1>, mp_size_t<1>, mp_quote<mp_list>>, mp_list<mp_list<X1>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2>, mp_size_t<1>, mp_quote<mp_list>>, mp_list<mp_list<X1>, mp_list<X2>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3>, mp_size_t<1>, mp_quote<mp_list>>, mp_list<mp_list<X1>, mp_list<X2>, mp_list<X3>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3, X4>, mp_size_t<1>, mp_quote<mp_list>>, mp_list<mp_list<X1>, mp_list<X2>, mp_list<X3>, mp_list<X4>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3, X4, X5>, mp_size_t<1>, mp_quote<mp_list>>, mp_list<mp_list<X1>, mp_list<X2>, mp_list<X3>, mp_list<X4>, mp_list<X5>>>));
    
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<>, mp_size_t<2>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1>, mp_size_t<2>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2>, mp_size_t<2>, mp_quote<mp_list>>, mp_list<mp_list<X1, X2>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3>, mp_size_t<2>, mp_quote<mp_list>>, mp_list<mp_list<X1, X2>, mp_list<X2, X3>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3, X4>, mp_size_t<2>, mp_quote<mp_list>>, mp_list<mp_list<X1, X2>, mp_list<X2, X3>, mp_list<X3, X4>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3, X4, X5>, mp_size_t<2>, mp_quote<mp_list>>, mp_list<mp_list<X1, X2>, mp_list<X2, X3>, mp_list<X3, X4>, mp_list<X4, X5>>>));
    
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<>, mp_size_t<3>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1>, mp_size_t<3>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2>, mp_size_t<3>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3>, mp_size_t<3>, mp_quote<mp_list>>, mp_list<mp_list<X1, X2, X3>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3, X4>, mp_size_t<3>, mp_quote<mp_list>>, mp_list<mp_list<X1, X2, X3>, mp_list<X2, X3, X4>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3, X4, X5>, mp_size_t<3>, mp_quote<mp_list>>, mp_list<mp_list<X1, X2, X3>, mp_list<X2, X3, X4>, mp_list<X3, X4, X5>>>));
    
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<>, mp_size_t<4>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1>, mp_size_t<4>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2>, mp_size_t<4>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3>, mp_size_t<4>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3, X4>, mp_size_t<4>, mp_quote<mp_list>>, mp_list<mp_list<X1, X2, X3, X4>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3, X4, X5>, mp_size_t<4>, mp_quote<mp_list>>, mp_list<mp_list<X1, X2, X3, X4>, mp_list<X2, X3, X4, X5>>>));
    
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<>, mp_size_t<5>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1>, mp_size_t<5>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2>, mp_size_t<5>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3>, mp_size_t<5>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3, X4>, mp_size_t<5>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3, X4, X5>, mp_size_t<5>, mp_quote<mp_list>>, mp_list<mp_list<X1, X2, X3, X4, X5>>>));
    
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<>, mp_size_t<6>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1>, mp_size_t<6>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2>, mp_size_t<6>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3>, mp_size_t<6>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3, X4>, mp_size_t<6>, mp_quote<mp_list>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3, X4, X5>, mp_size_t<6>, mp_quote<mp_list>>, mp_list<>>));
    
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<>, mp_size_t<2>, mp_quote<std::pair>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1>, mp_size_t<2>, mp_quote<std::pair>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2>, mp_size_t<2>, mp_quote<std::pair>>, mp_list<std::pair<X1, X2>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3>, mp_size_t<2>, mp_quote<std::pair>>, mp_list<std::pair<X1, X2>, std::pair<X2, X3>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3, X4>, mp_size_t<2>, mp_quote<std::pair>>, mp_list<std::pair<X1, X2>, std::pair<X2, X3>, std::pair<X3, X4>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list<X1, X2, X3, X4, X5>, mp_size_t<2>, mp_quote<std::pair>>, mp_list<std::pair<X1, X2>, std::pair<X2, X3>, std::pair<X3, X4>, std::pair<X4, X5>>>));
    
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<std::tuple<>, mp_size_t<2>, mp_quote<std::pair>>, std::tuple<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<std::tuple<X1>, mp_size_t<2>, mp_quote<std::pair>>, std::tuple<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<std::tuple<X1, X2>, mp_size_t<2>, mp_quote<std::pair>>, std::tuple<std::pair<X1, X2>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<std::tuple<X1, X2, X3>, mp_size_t<2>, mp_quote<std::pair>>, std::tuple<std::pair<X1, X2>, std::pair<X2, X3>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<std::tuple<X1, X2, X3, X4>, mp_size_t<2>, mp_quote<std::pair>>, std::tuple<std::pair<X1, X2>, std::pair<X2, X3>, std::pair<X3, X4>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<std::tuple<X1, X2, X3, X4, X5>, mp_size_t<2>, mp_quote<std::pair>>, std::tuple<std::pair<X1, X2>, std::pair<X2, X3>, std::pair<X3, X4>, std::pair<X4, X5>>>));
    
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list_c<int, 1, 2, 3, 4, 5, 6>, mp_size_t<1>, mp_quote<mp_plus>>, mp_list_c<int, 1, 2, 3, 4, 5, 6>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list_c<int, 1, 2, 3, 4, 5, 6>, mp_size_t<2>, mp_quote<mp_plus>>, mp_list_c<int, 3, 5, 7, 9, 11>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list_c<int, 1, 2, 3, 4, 5, 6>, mp_size_t<3>, mp_quote<mp_plus>>, mp_list_c<int, 6, 9, 12, 15>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list_c<int, 1, 2, 3, 4, 5, 6>, mp_size_t<4>, mp_quote<mp_plus>>, mp_list_c<int, 10, 14, 18>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list_c<int, 1, 2, 3, 4, 5, 6>, mp_size_t<5>, mp_quote<mp_plus>>, mp_list_c<int, 15, 20>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list_c<int, 1, 2, 3, 4, 5, 6>, mp_size_t<6>, mp_quote<mp_plus>>, mp_list_c<int, 21>>));
    
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list_c<int, 1, 2, 3, 4, 5, 6>, mp_size_t<1>, average>, mp_list_c<int, 1, 2, 3, 4, 5, 6>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list_c<int, 1, 2, 3, 4, 5, 6>, mp_size_t<2>, average>, mp_list_c<int, 1, 2, 3, 4, 5>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list_c<int, 1, 2, 3, 4, 5, 6>, mp_size_t<3>, average>, mp_list_c<int, 2, 3, 4, 5>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list_c<int, 1, 2, 3, 4, 5, 6>, mp_size_t<4>, average>, mp_list_c<int, 2, 3, 4>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list_c<int, 1, 2, 3, 4, 5, 6>, mp_size_t<5>, average>, mp_list_c<int, 3, 4>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_sliding_fold_q<mp_list_c<int, 1, 2, 3, 4, 5, 6>, mp_size_t<6>, average>, mp_list_c<int, 3>>));

    return boost::report_errors();
}
