/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// test_optional.cpp

// (C) Copyright 2004 Pavel Vozenilek
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// should pass compilation and execution

#include <cstddef> // NULL
#include <cstdio> // remove
#include <fstream>

#include <boost/config.hpp>

#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{
    using ::remove;
}
#endif

#include <boost/archive/archive_exception.hpp>
#include "test_tools.hpp"


struct A {
    int m_x;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* version */){
        ar & boost::serialization::make_nvp("x", m_x);
    };
    bool operator==(const A & rhs) const {
        return m_x == rhs.m_x;
    }
    // note that default constructor is not trivial
    A() :
        m_x(0)
    {}
    A(int x) :
        m_x(x)
    {}
};

// Optional is the class optional implementation you use
template<template<class> class Optional>
int test(){
    const char * testfile = boost::archive::tmpnam(NULL);
    BOOST_REQUIRE(NULL != testfile);

    const Optional<int> aoptional1;
    const Optional<int> aoptional2(123);
    const Optional<A> aoptional3;
    A a(1);
    const Optional<A> aoptional4(a);
    const Optional<A *> aoptional5;
    const Optional<A *> aoptional6(& a);
    {
        test_ostream os(testfile, TEST_STREAM_FLAGS);
        test_oarchive oa(os, TEST_ARCHIVE_FLAGS);
        oa << boost::serialization::make_nvp("aoptional1",aoptional1);
        oa << boost::serialization::make_nvp("aoptional2",aoptional2);
        oa << boost::serialization::make_nvp("aoptional3",aoptional3);
        oa << boost::serialization::make_nvp("aoptional4",aoptional4);
        oa << boost::serialization::make_nvp("aoptional5",aoptional5);
        oa << boost::serialization::make_nvp("aoptional6",aoptional6);
    }
    Optional<int> aoptional1a(999);
    Optional<int> aoptional2a;
    Optional<A> aoptional3a;
    Optional<A> aoptional4a;
    Optional<A *> aoptional5a;
    Optional<A *> aoptional6a;
    {
        test_istream is(testfile, TEST_STREAM_FLAGS);
        test_iarchive ia(is, TEST_ARCHIVE_FLAGS);
        ia >> boost::serialization::make_nvp("aoptional1",aoptional1a);
        ia >> boost::serialization::make_nvp("aoptional2",aoptional2a);
        ia >> boost::serialization::make_nvp("aoptional3",aoptional3a);
        ia >> boost::serialization::make_nvp("aoptional4",aoptional4a);
        ia >> boost::serialization::make_nvp("aoptional5",aoptional5a);
        ia >> boost::serialization::make_nvp("aoptional6",aoptional6a);
    }
    BOOST_CHECK(aoptional1 == aoptional1a);
    BOOST_CHECK(aoptional2 == aoptional2a);
    BOOST_CHECK(aoptional3 == aoptional3a);
    BOOST_CHECK(aoptional4 == aoptional4a);
    BOOST_CHECK(aoptional5 == aoptional5a);  // not initialized
    BOOST_CHECK(**aoptional6 == **aoptional6a);

    std::remove(testfile);
    return EXIT_SUCCESS;
}

#include <boost/serialization/optional.hpp>
#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
#include <optional>
#endif

int test_main( int /* argc */, char* /* argv */[] ){
    test<boost::optional>();
    #ifndef BOOST_NO_CXX17_HDR_OPTIONAL
    test<std::optional>();
    #endif
    return EXIT_SUCCESS;
}
