// Copyright 2017, 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#if defined(__clang__)
# pragma clang diagnostic ignored "-Wmismatched-tags"
#endif

#include <boost/container_hash/is_tuple_like.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <boost/config.hpp>
#include <boost/config/workaround.hpp>
#include <utility>

struct X
{
};

#if !defined(BOOST_NO_CXX11_HDR_TUPLE)

#include <tuple>

namespace user
{

struct Y
{
    int a;
    int b;
};

template<std::size_t I> int& get( Y& v );
template<std::size_t I> int const& get( Y const& v );

template<> int& get<0>( Y& v )
{
    return v.a;
}

template<> int const& get<0>( Y const& v )
{
    return v.a;
}

template<> int& get<1>( Y& v )
{
    return v.b;
}

template<> int const& get<1>( Y const& v )
{
    return v.b;
}

} // namespace user

namespace std
{

template<> struct tuple_size<user::Y>: std::integral_constant<std::size_t, 2>
{
};

} // namespace std

#endif // #if !defined(BOOST_NO_CXX11_HDR_TUPLE)

int main()
{
    using boost::container_hash::is_tuple_like;

    BOOST_TEST_TRAIT_FALSE((is_tuple_like<void>));
    BOOST_TEST_TRAIT_FALSE((is_tuple_like<void const>));

    BOOST_TEST_TRAIT_FALSE((is_tuple_like<int>));
    BOOST_TEST_TRAIT_FALSE((is_tuple_like<int const>));

    BOOST_TEST_TRAIT_FALSE((is_tuple_like<X>));
    BOOST_TEST_TRAIT_FALSE((is_tuple_like<X const>));

    BOOST_TEST_TRAIT_FALSE((is_tuple_like<int[2]>));
    BOOST_TEST_TRAIT_FALSE((is_tuple_like<int const [2]>));

#if !defined(BOOST_NO_CXX11_HDR_TUPLE) && !BOOST_WORKAROUND(BOOST_MSVC, <= 1800)

    BOOST_TEST_TRAIT_TRUE((is_tuple_like< std::pair<int, X> >));
    BOOST_TEST_TRAIT_TRUE((is_tuple_like< std::pair<int, X> const>));

    BOOST_TEST_TRAIT_TRUE((is_tuple_like< std::tuple<> >));
    BOOST_TEST_TRAIT_TRUE((is_tuple_like< std::tuple<> const>));

    BOOST_TEST_TRAIT_TRUE((is_tuple_like< std::tuple<X> >));
    BOOST_TEST_TRAIT_TRUE((is_tuple_like< std::tuple<X> const>));

    BOOST_TEST_TRAIT_TRUE((is_tuple_like< std::tuple<X, X> >));
    BOOST_TEST_TRAIT_TRUE((is_tuple_like< std::tuple<X, X> const>));

    BOOST_TEST_TRAIT_TRUE((is_tuple_like< std::tuple<X, X, X> >));
    BOOST_TEST_TRAIT_TRUE((is_tuple_like< std::tuple<X, X, X> const>));

    BOOST_TEST_TRAIT_TRUE((is_tuple_like<user::Y>));
    BOOST_TEST_TRAIT_TRUE((is_tuple_like<user::Y const>));

#endif

    return boost::report_errors();
}
