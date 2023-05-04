// Copyright 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/container_hash/is_described_class.hpp>
#include <boost/describe/class.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <string>

struct X1
{
};

struct X2
{
    int m;
};

struct Y1
{
};

BOOST_DESCRIBE_STRUCT( Y1, (), () )

struct Y2: Y1
{
    int m;
};

BOOST_DESCRIBE_STRUCT( Y2, (Y1), (m) )

int main()
{
    using boost::container_hash::is_described_class;

    BOOST_TEST_TRAIT_FALSE((is_described_class<void>));
    BOOST_TEST_TRAIT_FALSE((is_described_class<int>));
    BOOST_TEST_TRAIT_FALSE((is_described_class<X1>));
    BOOST_TEST_TRAIT_FALSE((is_described_class<X2>));
    BOOST_TEST_TRAIT_FALSE((is_described_class<int[2]>));
    BOOST_TEST_TRAIT_FALSE((is_described_class<std::string>));

#if defined(BOOST_DESCRIBE_CXX14)

    BOOST_TEST_TRAIT_TRUE((is_described_class<Y1>));
    BOOST_TEST_TRAIT_TRUE((is_described_class<Y2>));

#else

    BOOST_TEST_TRAIT_FALSE((is_described_class<Y1>));
    BOOST_TEST_TRAIT_FALSE((is_described_class<Y2>));

#endif

    return boost::report_errors();
}
