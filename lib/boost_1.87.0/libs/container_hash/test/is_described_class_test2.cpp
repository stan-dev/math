// Copyright 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/type_traits/integral_constant.hpp>
#include <boost/describe/class.hpp>

struct X
{
};

BOOST_DESCRIBE_STRUCT( X, (), () )

namespace boost
{
namespace container_hash
{

template<class T> struct is_described_class;
template<> struct is_described_class<X>: boost::false_type {};

} // namespace container_hash
} // namespace boost

#include <boost/container_hash/is_described_class.hpp>
#include <boost/core/lightweight_test_trait.hpp>

int main()
{
    using boost::container_hash::is_described_class;

    BOOST_TEST_TRAIT_FALSE((is_described_class<X>));

    return boost::report_errors();
}
