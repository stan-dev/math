// Copyright 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/container_hash/is_described_class.hpp>
#include <boost/describe/class.hpp>
#include <boost/core/lightweight_test_trait.hpp>

union Y1
{
};

BOOST_DESCRIBE_STRUCT( Y1, (), () )

union Y2
{
    int m1;
    float m2;
};

BOOST_DESCRIBE_STRUCT( Y2, (), (m1, m2) )

int main()
{
    using boost::container_hash::is_described_class;

    BOOST_TEST_TRAIT_FALSE((is_described_class<Y1>));
    BOOST_TEST_TRAIT_FALSE((is_described_class<Y2>));

    return boost::report_errors();
}
