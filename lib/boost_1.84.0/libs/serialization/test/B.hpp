#ifndef BOOST_SERIALIZATION_TEST_B_HPP
#define BOOST_SERIALIZATION_TEST_B_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// B.hpp

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com .
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <cstdlib> // for rand()
#include <cmath>
#include <boost/config.hpp>

#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{
    using ::rand;
}
#endif

#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/base_object.hpp>

#include "A.hpp"

///////////////////////////////////////////////////////
// Derived class test
class B : public A
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void save(Archive &ar, const unsigned int /* file_version */) const
    {
        // write any base class info to the archive
        ar << BOOST_SERIALIZATION_BASE_OBJECT_NVP(A);
    }

    template<class Archive>
    void load(Archive & ar, const unsigned int file_version)
    {
        // read any base class info to the archive
        ar >> BOOST_SERIALIZATION_BASE_OBJECT_NVP(A);
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()
    signed char s;
    unsigned char t;
    signed int u;
    unsigned int v;
    float w;
    double x;
public:
    B();
    virtual ~B(){};
    bool operator==(const B &rhs) const;
};

B::B()
{
}

BOOST_CLASS_VERSION(B, 2)

inline bool B::operator==(const B &rhs) const
{
    return A::operator==(rhs);
}

#endif // BOOST_SERIALIZATION_TEST_B_HPP
