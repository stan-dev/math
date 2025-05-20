// Copyright 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/type_traits/integral_constant.hpp>

struct X1
{
    char const* begin() const;
    char const* end() const;
};

struct X2
{
    char const* begin() const;
    char const* end() const;
};

namespace boost
{
namespace container_hash
{

template<class T> struct is_range;
template<> struct is_range<X2>: boost::false_type {};

} // namespace container_hash
} // namespace boost

#include <boost/container_hash/is_range.hpp>
#include <boost/core/lightweight_test_trait.hpp>

int main()
{
    using boost::container_hash::is_range;

#if !defined(BOOST_NO_CXX11_DECLTYPE) && !defined(BOOST_NO_SFINAE_EXPR) && !BOOST_WORKAROUND(BOOST_GCC, < 40700)

    BOOST_TEST_TRAIT_TRUE((is_range<X1>));

#endif

    BOOST_TEST_TRAIT_FALSE((is_range<X2>));

    return boost::report_errors();
}
