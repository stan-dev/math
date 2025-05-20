// Copyright 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/type_traits/integral_constant.hpp>
#include <cstddef>

struct X1
{
    char const* begin() const;
    char const* end() const;
    typedef void hasher;
};

struct X2
{
    char const* begin() const;
    char const* end() const;
    typedef void hasher;
};

struct X3
{
    char const* begin() const;
    char const* end() const;
    typedef void hasher;
};

namespace boost
{
namespace container_hash
{

template<class T> struct is_unordered_range;
template<> struct is_unordered_range<X2>: boost::false_type {};

template<class T> struct is_range;
template<> struct is_range<X3>: boost::false_type {};

} // namespace container_hash
} // namespace boost

#include <boost/container_hash/is_unordered_range.hpp>
#include <boost/core/lightweight_test_trait.hpp>

int main()
{
    using boost::container_hash::is_unordered_range;

#if !defined(BOOST_NO_CXX11_DECLTYPE) && !defined(BOOST_NO_SFINAE_EXPR) && !BOOST_WORKAROUND(BOOST_GCC, < 40700)

    BOOST_TEST_TRAIT_TRUE((is_unordered_range<X1>));

#endif

    BOOST_TEST_TRAIT_FALSE((is_unordered_range<X2>));
    BOOST_TEST_TRAIT_FALSE((is_unordered_range<X3>));

    return boost::report_errors();
}
