// Copyright 2015, 2023 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/integral.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <type_traits>

int main()
{
    using boost::mp11::mp_list;
    using boost::mp11::mp_iota;
    using boost::mp11::mp_iota_c;
    using boost::mp11::mp_int;
    using boost::mp11::mp_size_t;

    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_iota_c<0, 7>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_iota_c<1, 7>, mp_list<mp_size_t<7>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_iota_c<2, 7>, mp_list<mp_size_t<7>, mp_size_t<8>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_iota_c<3, 7>, mp_list<mp_size_t<7>, mp_size_t<8>, mp_size_t<9>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_iota_c<4, 7>, mp_list<mp_size_t<7>, mp_size_t<8>, mp_size_t<9>, mp_size_t<10>>>));

    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_iota<mp_size_t<0>, mp_int<4>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_iota<mp_size_t<1>, mp_int<4>>, mp_list<mp_size_t<4>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_iota<mp_size_t<2>, mp_int<4>>, mp_list<mp_size_t<4>, mp_size_t<5>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_iota<mp_size_t<3>, mp_int<4>>, mp_list<mp_size_t<4>, mp_size_t<5>, mp_size_t<6>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_iota<mp_size_t<4>, mp_int<4>>, mp_list<mp_size_t<4>, mp_size_t<5>, mp_size_t<6>, mp_size_t<7>>>));

    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_iota<mp_int<0>, mp_int<-2>>, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_iota<mp_int<1>, mp_int<-2>>, mp_list<mp_int<-2>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_iota<mp_int<2>, mp_int<-2>>, mp_list<mp_int<-2>, mp_int<-1>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_iota<mp_int<3>, mp_int<-2>>, mp_list<mp_int<-2>, mp_int<-1>, mp_int<0>>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_iota<mp_int<4>, mp_int<-2>>, mp_list<mp_int<-2>, mp_int<-1>, mp_int<0>, mp_int<+1>>>));

    return boost::report_errors();
}
