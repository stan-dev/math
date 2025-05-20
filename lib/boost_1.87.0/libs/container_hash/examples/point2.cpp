
// Copyright 2005 Daniel James.
// Copyright 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Force use of assert.
#if defined(NDEBUG)
#undef NDEBUG
#endif

#include <boost/container_hash/hash.hpp>
#include <boost/describe/class.hpp>
#include <boost/describe/operators.hpp>
#include <cassert>

// This example illustrates how to use Boost.Describe to obtain
// automatic boost::hash support. For full details see the hash
// tutorial.

#if defined(BOOST_DESCRIBE_CXX14)

class point
{
    int x;
    int y;

    BOOST_DESCRIBE_CLASS(point, (), (), (), (x, y))

public:

    point() : x(0), y(0) {}
    point(int x, int y) : x(x), y(y) {}
};

using boost::describe::operators::operator==;
using boost::describe::operators::operator!=;

int main()
{
    boost::hash<point> point_hasher;

    point p1(0, 0);
    point p2(1, 2);
    point p3(4, 1);
    point p4 = p1;

    assert(point_hasher(p1) == point_hasher(p4));

    // These tests could legally fail, but if they did it'd be a pretty bad
    // hash function.
    assert(point_hasher(p1) != point_hasher(p2));
    assert(point_hasher(p1) != point_hasher(p3));
}

#else

#include <cstdio>

int main()
{
    std::puts( "This example requires C++14." );
}

#endif
