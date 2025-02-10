// Copyright 2024 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/container_hash/is_range.hpp>
#include <boost/core/lightweight_test_trait.hpp>

struct X1
{
    void* begin() const;
    void* end() const;
};

struct X2
{
    void const* begin() const;
    void const* end() const;
};

int main()
{
    using boost::container_hash::is_range;

    BOOST_TEST_TRAIT_FALSE((is_range<X1>));
    BOOST_TEST_TRAIT_FALSE((is_range<X2>));

    return boost::report_errors();
}
