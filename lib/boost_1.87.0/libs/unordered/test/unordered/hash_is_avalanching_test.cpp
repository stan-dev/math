// Copyright 2022 Peter Dimov
// Copyright 2024 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// an imitation of a third-party header specializing hash_is_avalanching
// (boost/container_hash/hash.hpp is an example doing that)

#include <boost/type_traits/integral_constant.hpp>

struct X3
{
};

namespace boost
{
namespace unordered
{

    template<class T> struct hash_is_avalanching;
    template<> struct hash_is_avalanching< ::X3 >: boost::true_type {};

} // namespace unordered
} // namespace boost

//

#include <boost/unordered/hash_traits.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <type_traits>

struct X1
{
};

struct X2
{
    typedef void is_avalanching;
};

struct X4
{
    using is_avalanching = std::false_type;
};

struct X5
{
    using is_avalanching = std::true_type;
};

struct X6
{
    using is_avalanching = boost::false_type;
};

struct X7
{
    using is_avalanching = boost::true_type;
};


int main()
{
    using boost::unordered::hash_is_avalanching;

    BOOST_TEST_TRAIT_FALSE((hash_is_avalanching<X1>));
    BOOST_TEST_TRAIT_TRUE((hash_is_avalanching<X2>));
    BOOST_TEST_TRAIT_TRUE((hash_is_avalanching<X3>));
    BOOST_TEST_TRAIT_FALSE((hash_is_avalanching<X4>));
    BOOST_TEST_TRAIT_TRUE((hash_is_avalanching<X5>));
    BOOST_TEST_TRAIT_FALSE((hash_is_avalanching<X6>));
    BOOST_TEST_TRAIT_TRUE((hash_is_avalanching<X7>));

    return boost::report_errors();
}
