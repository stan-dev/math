// Copyright 2005-2009 Daniel James.
// Copyright 2021 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/container_hash/hash.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/core/type_name.hpp>
#include <set>

#if defined(BOOST_MSVC)
# pragma warning(disable: 4244) // conversion from int to float
#endif

#if defined(__GNUC__) || defined(__clang__)
# pragma GCC diagnostic ignored "-Wconversion"
# pragma GCC diagnostic ignored "-Wsign-conversion"
#endif

// This test checks for collisions in a small range of numbers

template<class T, int M>
void collision_test_()
{
    std::set<std::size_t> hashes;

    for( int i = -128; i <= 127; ++i )
    {
        hashes.insert( boost::hash<T>()( i * M ) );
    }

    BOOST_TEST_EQ( hashes.size(), 256u );
}

template <class T>
void collision_test()
{
    collision_test_<T, 1>();
    collision_test_<T, 2>();
    collision_test_<T, 3>();
    collision_test_<T, 4>();
    collision_test_<T, 5>();
    collision_test_<T, 8>();
    collision_test_<T, 10>();
    collision_test_<T, 16>();
    collision_test_<T, 32>();
    collision_test_<T, 64>();
    collision_test_<T, 100>();
    collision_test_<T, 128>();
}

#define TEST(type) std::cerr << "Testing: " #type " (" << boost::core::type_name<type>() << ")\n"; collision_test<type>();

int main()
{
#ifndef BOOST_NO_INTRINSIC_WCHAR_T
    TEST(wchar_t)
#endif
#ifndef BOOST_NO_CXX11_CHAR16_T
    TEST(char16_t)
#endif
#ifndef BOOST_NO_CXX11_CHAR32_T
    TEST(char32_t)
#endif
    TEST(short)
    TEST(unsigned short)
    TEST(int)
    TEST(unsigned int)
    TEST(long)
    TEST(unsigned long)

#if !defined(BOOST_NO_LONG_LONG)
    TEST(boost::long_long_type)
    TEST(boost::ulong_long_type)
#endif

#if defined(BOOST_HAS_INT128)
    TEST(boost::int128_type)
    TEST(boost::uint128_type)
#endif

    TEST(float)
    TEST(double)
    TEST(long double)

    TEST(std::size_t)
    TEST(std::ptrdiff_t)

    return boost::report_errors();
}
