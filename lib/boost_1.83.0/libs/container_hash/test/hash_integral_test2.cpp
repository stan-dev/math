// Copyright 2005-2009 Daniel James.
// Copyright 2021 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/container_hash/hash.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/core/type_name.hpp>
#include <boost/type_traits/is_signed.hpp>
#include <boost/type_traits/make_unsigned.hpp>
#include <boost/static_assert.hpp>
#include <boost/config.hpp>

// This test checks that values representable in a signed
// and the corresponding unsigned type hash to the same value

template<class T>
void signed_unsigned_test()
{
    BOOST_STATIC_ASSERT( boost::is_signed<T>::value );

    typedef typename boost::make_unsigned<T>::type U;

    T x = std::numeric_limits<T>::max();

    do
    {
        BOOST_TEST_EQ( boost::hash<T>()( x ), boost::hash<U>()( static_cast<U>( x ) ) );
        x /= 3;
    }
    while( x > 0 );
}

#define TEST(type) std::cerr << "Testing: " #type " (" << boost::core::type_name<type>() << ")\n"; signed_unsigned_test<type>();

int main()
{
    TEST(signed char)
    TEST(short)
    TEST(int)
    TEST(long)

#if !defined(BOOST_NO_LONG_LONG)
    TEST(boost::long_long_type)
#endif

#if defined(BOOST_HAS_INT128)
    TEST(boost::int128_type)
#endif

    return boost::report_errors();
}
