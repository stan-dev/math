// Copyright 2024 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/container_hash/hash.hpp>
#include <boost/core/lightweight_test.hpp>

struct X1
{
    void* begin() const { return nullptr; }
    void* end() const { return nullptr; }
};

std::size_t hash_value( X1 const& )
{
    return 1;
}

struct X2
{
    void const* begin() const { return nullptr; }
    void const* end() const { return nullptr; }
};

std::size_t hash_value( X2 const& )
{
    return 2;
}

int main()
{
    X1 x1;
    BOOST_TEST_EQ( boost::hash<X1>()( x1 ), 1u );

    X2 x2;
    BOOST_TEST_EQ( boost::hash<X2>()( x2 ), 2u );

    return boost::report_errors();
}
