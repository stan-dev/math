// Copyright 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/container_hash/is_contiguous_range.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <cstddef>

struct X1
{
    char const* begin() const;
    char const* end() const;
    std::size_t size() const;

    char data[16];
};

int main()
{
    using boost::container_hash::is_contiguous_range;

    BOOST_TEST_TRAIT_FALSE((is_contiguous_range<X1>));

    return boost::report_errors();
}
