///////////////////////////////////////////////////////////////
//  Copyright 2024 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt

#include <iostream>
#include <boost/multiprecision/cpp_int.hpp>
#include "test.hpp"

template <class T>
void test_subtract_underflow()
{
    T val = 42;

    int sub = 43;

    BOOST_CHECK_THROW(static_cast<T>(val - sub), std::range_error);
    BOOST_CHECK_THROW(static_cast<T>(val - static_cast<T>(sub)), std::range_error);
    BOOST_CHECK_THROW(static_cast<T>(val - static_cast<std::uintmax_t>(sub)), std::range_error);
    BOOST_CHECK_THROW(static_cast<T>(val - static_cast<std::intmax_t>(sub)), std::range_error);

    BOOST_IF_CONSTEXPR((std::numeric_limits<T>::digits < 500) && std::numeric_limits<T>::digits)
    {
       val <<= std::numeric_limits<T>::digits - 10;
    }
    else
       val <<= 500;

    T sub2 = val + 1;

    BOOST_CHECK_THROW(static_cast<T>(val - sub2), std::range_error);
}

template <unsigned int N>
void test_overflow() {
   using uint_t = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<N, N, boost::multiprecision::unsigned_magnitude, boost::multiprecision::cpp_int_check_type::checked, void>>;

   uint_t value = (std::numeric_limits<uint_t>::max)();
   BOOST_CHECK_THROW(value += 1, std::overflow_error);
   //
   // We don't care what the value is, but it must be sane and normalized
   // to be within the range of the type under test:
   //
   BOOST_CHECK(value <= (std::numeric_limits<uint_t>::max)());
}

template <unsigned int N>
struct TestOverflow {
   static void run() {
      test_overflow<N>();
      TestOverflow<N + 8>::run();
   }
};

template <>
struct TestOverflow<1032> {
   static void run() {}
};



int main()
{
    using checked_uint64_t = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<64, 64, boost::multiprecision::unsigned_magnitude, boost::multiprecision::checked>>;

    test_subtract_underflow<checked_uint64_t>();
    test_subtract_underflow<boost::multiprecision::checked_uint128_t>();
    test_subtract_underflow<boost::multiprecision::checked_uint256_t>();
    test_subtract_underflow<boost::multiprecision::checked_uint1024_t>();

   TestOverflow<8>::run();

    return boost::report_errors();
}

