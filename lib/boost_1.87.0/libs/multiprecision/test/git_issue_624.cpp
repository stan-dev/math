///////////////////////////////////////////////////////////////
//  Copyright 2024 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt

#include <iostream>
#include <boost/multiprecision/cpp_int.hpp>
#include "test.hpp"

template <class T>
void test_neg_divide_by_zero(std::true_type const&)
{
    T val = -42;
    BOOST_CHECK_THROW(static_cast<T>(val % 0), std::overflow_error);
    BOOST_CHECK_THROW(static_cast<T>(val % static_cast<T>(0)), std::overflow_error);
    BOOST_CHECK_THROW(static_cast<T>(val % static_cast<std::uintmax_t>(0)), std::overflow_error);
    BOOST_CHECK_THROW(static_cast<T>(val % static_cast<std::intmax_t>(0)), std::overflow_error);
}
template <class T>
void test_neg_divide_by_zero(std::false_type const&)
{
}

template <class T>
void test_divide_by_zero()
{
    T val = 42;

    BOOST_CHECK_THROW(static_cast<T>(val % 0), std::overflow_error);
    BOOST_CHECK_THROW(static_cast<T>(val % static_cast<T>(0)), std::overflow_error);
    BOOST_CHECK_THROW(static_cast<T>(val % static_cast<std::uintmax_t>(0)), std::overflow_error);
    BOOST_CHECK_THROW(static_cast<T>(val % static_cast<std::intmax_t>(0)), std::overflow_error);

    BOOST_IF_CONSTEXPR((std::numeric_limits<T>::digits < 500) && std::numeric_limits<T>::digits)
    {
       val <<= std::numeric_limits<T>::digits - 10;
    }
    else
       val <<= 500;

    BOOST_CHECK_THROW(static_cast<T>(val % 0), std::overflow_error);
    BOOST_CHECK_THROW(static_cast<T>(val % static_cast<T>(0)), std::overflow_error);
    BOOST_CHECK_THROW(static_cast<T>(val % static_cast<std::uintmax_t>(0)), std::overflow_error);
    BOOST_CHECK_THROW(static_cast<T>(val % static_cast<std::intmax_t>(0)), std::overflow_error);

    using tag_type = std::integral_constant<bool, std::numeric_limits<T>::is_signed>;

    test_neg_divide_by_zero<T>(tag_type());
}

int main()
{
    using int64_t = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<64, 64, boost::multiprecision::signed_magnitude>>;
    using uint64_t = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<64, 64, boost::multiprecision::unsigned_magnitude>>;
    using checked_int64_t = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<64, 64, boost::multiprecision::signed_magnitude, boost::multiprecision::checked>>;
    using checked_uint64_t = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<64, 64, boost::multiprecision::unsigned_magnitude, boost::multiprecision::checked>>;

    test_divide_by_zero<int64_t>();
    test_divide_by_zero<boost::multiprecision::int128_t>();
    test_divide_by_zero<boost::multiprecision::int256_t>();
    test_divide_by_zero<boost::multiprecision::int1024_t>();

    test_divide_by_zero<uint64_t>();
    test_divide_by_zero<boost::multiprecision::uint128_t>();
    test_divide_by_zero<boost::multiprecision::uint256_t>();
    test_divide_by_zero<boost::multiprecision::uint1024_t>();

    test_divide_by_zero<checked_int64_t>();
    test_divide_by_zero<boost::multiprecision::checked_int128_t>();
    test_divide_by_zero<boost::multiprecision::checked_int256_t>();
    test_divide_by_zero<boost::multiprecision::checked_int1024_t>();

    test_divide_by_zero<uint64_t>();
    test_divide_by_zero<boost::multiprecision::checked_uint128_t>();
    test_divide_by_zero<boost::multiprecision::checked_uint256_t>();
    test_divide_by_zero<boost::multiprecision::checked_uint1024_t>();

    test_divide_by_zero<boost::multiprecision::cpp_int>();

    return boost::report_errors();
}

