// Copyright 2024 Matt Borland
// Copyright 2024 Junekey Jeon
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/cppalliance/charconv/issues/158

#include <boost/charconv.hpp>
#include <boost/core/lightweight_test.hpp>

void test_values_with_negative_exp()
{
    char buffer[256];

    double d = 1e-13;
    std::memset(buffer, '\0', sizeof(buffer));
    auto res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::scientific, 50);
    *res.ptr = '\0';

    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1.00000000000000003037374556340037091360347168422784e-13");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
                         boost::charconv::chars_format::fixed, 50);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "0.00000000000010000000000000000303737455634003709136");

    d = 1e-15;
    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
                         boost::charconv::chars_format::scientific, 50);
    *res.ptr = '\0';

    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1.00000000000000007770539987666107923830718560119502e-15");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::fixed, 50);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "0.00000000000000100000000000000007770539987666107924");

    d = 1e-17;

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::scientific, 50);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1.00000000000000007154242405462192450852805618492325e-17");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::fixed, 50);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "0.00000000000000001000000000000000071542424054621925");

    d = 1e-19;

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), d,
                         boost::charconv::chars_format::scientific, 50);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "9.99999999999999975245926835260131855729159055676882e-20");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), d,
                         boost::charconv::chars_format::fixed, 50);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "0.00000000000000000009999999999999999752459268352601");

    d = 1e-21;

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), d,
                         boost::charconv::chars_format::scientific, 50);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "9.99999999999999907537452227896371396729934511675531e-22");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), d,
                         boost::charconv::chars_format::fixed, 50);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "0.00000000000000000000099999999999999990753745222790");
}

#ifndef BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE
void test_long_double_with_negative_exp()
{
    char buffer[256];
    long double d = 1e-15L;
    auto res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), d,
                         boost::charconv::chars_format::scientific, 50);
    *res.ptr = '\0';

    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "9.99999999999999999994126620636112567498475617893198e-16");

    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), d,
                         boost::charconv::chars_format::fixed, 50);
    *res.ptr = '\0';
    BOOST_TEST(res);
    // BOOST_TEST_CSTR_EQ(buffer, "0.00000000000000099999999999999999999412662063611257");
    BOOST_TEST_CSTR_EQ(buffer, "0.00000000000000100000000000000000000000000000000000");

    d = 1e-17L;

    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), d,
                         boost::charconv::chars_format::scientific, 50);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "9.99999999999999999997135886174217623518875583428487e-18");

    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), d,
                         boost::charconv::chars_format::fixed, 50);
    *res.ptr = '\0';
    BOOST_TEST(res);
    // BOOST_TEST_CSTR_EQ(buffer, "0.00000000000000000999999999999999999997135886174218");
    BOOST_TEST_CSTR_EQ(buffer, "0.00000000000000001000000000000000000000000000000000");
}
#endif

void test_values_with_positive_exp()
{
    char buffer[256];
    double d = 1e15;
    auto res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::scientific, 50);
    *res.ptr = '\0';

    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1.00000000000000000000000000000000000000000000000000e+15");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::fixed, 50);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1000000000000000.00000000000000000000000000000000000000000000000000");

    d = 1e17;

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::scientific, 50);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1.00000000000000000000000000000000000000000000000000e+17");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::fixed, 50);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "100000000000000000.00000000000000000000000000000000000000000000000000");
}

void test_round_9()
{
    char buffer[256];
    double d = 999999999.999999;
    auto res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::scientific, 10);
    *res.ptr = '\0';

    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1.0000000000e+09");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::fixed, 3);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1000000000.000");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::general, 9);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1e+09");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::general, 10);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1000000000");

    d = 999999.999999;

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::scientific, 10);
    *res.ptr = '\0';

    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1.0000000000e+06");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::fixed, 3);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1000000.000");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::general, 6);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1e+06");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::general, 7);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1000000");

    d = 9.999999;

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::scientific, 5);
    *res.ptr = '\0';

    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1.00000e+01");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::fixed, 3);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "10.000");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::general, 6);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "10");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::general, 7);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "9.999999");

    d = 0.9999999;

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::scientific, 5);
    *res.ptr = '\0';

    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1.00000e+00");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::fixed, 3);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1.000");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::general, 6);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::general, 7);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "0.9999999");

    d = 0.0009999999;

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::scientific, 5);
    *res.ptr = '\0';

    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1.00000e-03");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::fixed, 3);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "0.001");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::general, 6);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "0.001");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::general, 7);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "0.0009999999");

    d = 0.00009999999;

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::scientific, 5);
    *res.ptr = '\0';

    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1.00000e-04");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::fixed, 3);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "0.000");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::general, 6);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "0.0001");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::general, 7);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "9.999999e-05");

    d = 0.00000009999999;

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::scientific, 5);
    *res.ptr = '\0';

    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1.00000e-07");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::fixed, 3);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "0.000");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::general, 6);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1e-07");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::general, 7);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "9.999999e-08");
}

void test_zero()
{
    char buffer[256];
    double d = 0;
    auto res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::scientific, 50);
    *res.ptr = '\0';

    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "0.00000000000000000000000000000000000000000000000000e+00");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::fixed, 50);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "0.00000000000000000000000000000000000000000000000000");

    std::memset(buffer, '\0', sizeof(buffer));
    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, d,
        boost::charconv::chars_format::general, 50);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "0");
}

#ifndef BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE
void test_long_double_with_positive_exp()
{
    char buffer[256];
    long double d = 1e15L;
    auto res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), d,
                         boost::charconv::chars_format::scientific, 50);
    *res.ptr = '\0';

    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1.00000000000000000000000000000000000000000000000000e+15");

    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), d,
                         boost::charconv::chars_format::fixed, 50);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1000000000000000.00000000000000000000000000000000000000000000000000");

    d = 1e17L;

    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), d,
                         boost::charconv::chars_format::scientific, 50);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "1.00000000000000000000000000000000000000000000000000e+17");

    res = boost::charconv::to_chars(buffer, buffer + sizeof(buffer), d,
                         boost::charconv::chars_format::fixed, 50);
    *res.ptr = '\0';
    BOOST_TEST(res);
    BOOST_TEST_CSTR_EQ(buffer, "100000000000000000.00000000000000000000000000000000000000000000000000");
}
#endif

template <typename T>
void test_spot_value(T value, int precision, const char* result, boost::charconv::chars_format fmt = boost::charconv::chars_format::fixed)
{
    char buffer[256];
    auto r = boost::charconv::to_chars(buffer, buffer + sizeof(buffer) - 1, value, fmt, precision);
    *r.ptr = '\0';

    BOOST_TEST(r);
    BOOST_TEST_CSTR_EQ(buffer, result);
}

int main()
{
    test_values_with_negative_exp();
    test_values_with_positive_exp();
    test_zero();
    test_round_9();

    #if BOOST_CHARCONV_LDBL_BITS == 80
    test_long_double_with_negative_exp();
    test_long_double_with_positive_exp();
    #endif

    // Found during random testing in to_chars_float_STL_comp
    //test_spot_value(27057.375F, 49, "27057.3750000000000000000000000000000000000000000000000");
    //test_spot_value(-38347.10547F, 49, "-38347.1054687500000000000000000000000000000000000000000");
    //test_spot_value(12043.7270408630284, 49, "12043.727040863028378225862979888916015625", boost::charconv::chars_format::general);

    #if BOOST_CHARCONV_LDBL_BITS == 80
    test_spot_value(-37512.44400929347152385L, 17, "-3.75124440092934715e+04", boost::charconv::chars_format::scientific);
    test_spot_value(69260.5792027112860012L, 17, "6.92605792027112860e+04", boost::charconv::chars_format::scientific);
    test_spot_value(7.420390267538564899041L,  13, "7.4203902675386e+00", boost::charconv::chars_format::scientific);
    test_spot_value(-72911.41743909836600324L, 15, "-7.291141743909837e+04", boost::charconv::chars_format::scientific);
    #endif

    return boost::report_errors();
}
