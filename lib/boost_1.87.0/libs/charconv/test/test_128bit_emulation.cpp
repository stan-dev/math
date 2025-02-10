// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/config.hpp>

#ifdef BOOST_HAS_INT128

// We need to define these operator<< overloads before
// including boost/core/lightweight_test.hpp, or they
// won't be visible to BOOST_TEST_EQ
// LCOV_EXCL_START

#include <ostream>

static char* mini_to_chars( char (&buffer)[ 64 ], boost::uint128_type v )
{
    char* p = buffer + 64;
    *--p = '\0';

    do
    {
        *--p = "0123456789"[ v % 10 ];
        v /= 10;
    }
    while ( v != 0 );

    return p;
}

std::ostream& operator<<( std::ostream& os, boost::uint128_type v )
{
    char buffer[ 64 ];

    os << mini_to_chars( buffer, v );
    return os;
}

std::ostream& operator<<( std::ostream& os, boost::int128_type v )
{
    char buffer[ 64 ];
    char* p;

    if( v >= 0 )
    {
        p = mini_to_chars( buffer, static_cast<boost::uint128_type>(v) );
    }
    else
    {
        p = mini_to_chars( buffer, -static_cast<boost::uint128_type>(v) );
        *--p = '-';
    }

    os << p;
    return os;
}

// LCOV_EXCL_STOP

#endif // #ifdef BOOST_HAS_INT128

#include <boost/charconv/detail/emulated128.hpp>
#include <boost/core/lightweight_test.hpp>
#include <limits>
#include <iostream>
#include <climits>
#include <cstdint>

using boost::charconv::detail::uint128;
using boost::charconv::detail::trivial_uint128;

template <typename T>
void test_relational_operators(T val = (std::numeric_limits<T>::max)())
{
    uint128 test_val = UINT64_MAX;
    test_val += 1;

    BOOST_TEST(test_val > val);
    BOOST_TEST(!(test_val < val));
    BOOST_TEST(!(test_val == val));
    BOOST_TEST(test_val != val);

    const uint128 equal_val = val;

    BOOST_TEST(!(equal_val > val));
    BOOST_TEST(equal_val >= val);
    BOOST_TEST(!(equal_val < val));
    BOOST_TEST(equal_val <= val);
    BOOST_TEST(equal_val == val);
    BOOST_TEST(!(equal_val != val));

    int negative_val = -100;

    BOOST_TEST(test_val > negative_val);
    BOOST_TEST(!(test_val < negative_val));
    BOOST_TEST(!(test_val == negative_val));
    BOOST_TEST(test_val != negative_val);
}

void test_arithmetic_operators()
{
    // Only using low word
    const auto fixed_val = UINT64_MAX / 2;
    uint128 test_val = fixed_val;
    BOOST_TEST(test_val / 2 == UINT64_MAX / 4);
    BOOST_TEST(test_val + 1 == fixed_val + 1);
    test_val++;
    BOOST_TEST(test_val == fixed_val + 1);
    BOOST_TEST(test_val % fixed_val == 1);
    test_val--;
    BOOST_TEST(test_val == fixed_val);
    BOOST_TEST(test_val % fixed_val == 0);
    BOOST_TEST(test_val / fixed_val == 1);


    test_val = 2;
    std::uint64_t comp_val = 1;
	while (test_val < UINT64_MAX)
    {
        comp_val *= 2;
		if(!BOOST_TEST(test_val == comp_val))
		{
            // LCOV_EXCL_START
            std::cerr << "Target: " << comp_val
                      << "\ntest_val: " << test_val.low << std::endl;
            // LCOV_EXCL_STOP
		}
        test_val *= 2;
    }

    // And back down
    while (test_val >= 2)
    {
        test_val /= 2;
        if(!BOOST_TEST(test_val == comp_val))
        {
            // LCOV_EXCL_START
            std::cerr << "Target: " << comp_val
                      << "\ntest_val: " << test_val.low << std::endl;
            // LCOV_EXCL_STOP
        }
        comp_val /= 2;
    }


    // Add the high word
    uint128 test_high_word = UINT64_MAX;
    ++test_high_word;
    BOOST_TEST(test_high_word.high == 1 && test_high_word.low == 0);
    --test_high_word;

	#ifdef BOOST_CHARCONV_HAS_INT128
    boost::uint128_type reference = UINT64_MAX;
    BOOST_TEST(test_high_word == reference);

    for (int i = 0; i < 63; ++i)
    {
        if(!BOOST_TEST(test_high_word == reference))
        {
            // LCOV_EXCL_START
            std::cerr << "i: " << i
                      << "\nTarget: " << reference
                      << "\ntest_val: " << test_high_word.high << " " << test_high_word.low << std::endl;
            // // LCOV_EXCL_STOP
        }
        test_high_word *= 2;
        reference *= 2;
    }

    while (test_high_word >= 2)
    {
        BOOST_TEST(test_high_word == reference);
        test_high_word /= 2;
        reference /= 2;
    }

	#endif
}

void test_bitwise_operators()
{
    #ifdef BOOST_CHARCONV_HAS_INT128
    boost::uint128_type ref = UINT64_MAX;
    uint128 test_val = UINT64_MAX;

    ref <<= 1;
    test_val <<= 1;
    BOOST_TEST(test_val == ref);

    ref >>= 2;
    test_val >>= 2;
    BOOST_TEST(test_val == ref);

    BOOST_TEST((test_val | 1) == (ref | 1));
    BOOST_TEST((test_val & 1) == (ref & 1));
    BOOST_TEST(~test_val == ~ref);
    #endif
}

void test_memcpy()
{
    #if defined(BOOST_CHARCONV_HAS_FLOAT128) && defined(BOOST_CHARCONV_HAS_INT128)
    __float128 fval = 1e4000Q;
    boost::uint128_type ref;
    trivial_uint128 cpyval;

    std::memcpy(&cpyval, &fval, sizeof(fval));
    std::memcpy(&ref, &fval, sizeof(fval));

    uint128 test_val = cpyval;

    BOOST_TEST(test_val == ref);
    #endif
}

void test_limits()
{
    BOOST_TEST(std::numeric_limits<uint128>::min() == 0);
    BOOST_TEST(std::numeric_limits<uint128>::max() == uint128(UINT64_MAX, UINT64_MAX));
    BOOST_TEST(std::numeric_limits<uint128>::is_signed == false);
    BOOST_TEST(std::numeric_limits<uint128>::is_integer == true);
    BOOST_TEST(std::numeric_limits<uint128>::digits == CHAR_BIT * sizeof(std::uint64_t) * 2);

    // Max value is 340,282,366,920,938,463,463,374,607,431,768,211,455 (39 digits) so 38 without change
    BOOST_TEST(std::numeric_limits<uint128>::digits10 == 38);
}

int main()
{
    test_relational_operators<char>();
    test_relational_operators<signed char>();
    test_relational_operators<short>();
    test_relational_operators<int>();
    test_relational_operators<long>();
    test_relational_operators<long long>();
    test_relational_operators<unsigned char>();
    test_relational_operators<unsigned short>();
    test_relational_operators<unsigned>();
    test_relational_operators<unsigned long>();
    test_relational_operators<unsigned long long>();

    test_arithmetic_operators();

    test_bitwise_operators();

    test_memcpy();

    test_limits();

    return boost::report_errors();
}
