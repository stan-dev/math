//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
// Copyright (c) 2021-2023 Alexander Grund
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include "../src/boost/locale/shared/mo_lambda.hpp"
#include "boostLocale/test/unit_test.hpp"
#include <iostream>
#include <limits>
#include <random>

template<typename T>
T getRandValue(const T min, const T max)
{
    static std::mt19937 gen{std::random_device{}()};
    std::uniform_int_distribution<T> distrib(min, max);
    return distrib(gen);
}

template<typename T>
void test_plural_expr_rand(const T& ref, const char* expr)
{
    constexpr auto minVal = std::numeric_limits<long long>::min() / 1024;
    constexpr auto maxVal = std::numeric_limits<long long>::max() / 1024;
    const auto ptr = boost::locale::gnu_gettext::lambda::compile(expr);
    TEST(ptr);
    constexpr int number_of_tries = 256;
    for(int i = 0; i < number_of_tries; i++) {
        const auto n = getRandValue(minVal, maxVal);
        const auto result = ptr(n);
        const auto refResult = ref(n);
        if(result != refResult) {
            std::cerr << "Expression: " << expr << "; n=" << n << '\n'; // LCOV_EXCL_LINE
            TEST_EQ(result, refResult);                                 // LCOV_EXCL_LINE
        }
    }
}

void test_plural_expr()
{
    using boost::locale::gnu_gettext::lambda::compile;
#define COMPILE_PLURAL_EXPR(expr) \
    []() {                        \
        auto ptr = compile(expr); \
        TEST(ptr);                \
        return ptr;               \
    }()
#define TEST_EQ_EXPR(expr, rhs) test_eq_impl(COMPILE_PLURAL_EXPR(expr)(0), rhs, expr, __LINE__)
    // Number only
    TEST_EQ_EXPR("0", 0);
    TEST_EQ_EXPR("42", 42);
    BOOST_LOCALE_START_CONST_CONDITION
    if(sizeof(long) >= 4) {
        BOOST_LOCALE_END_CONST_CONDITION
        TEST_EQ_EXPR("2147483647", 2147483647); // largest signed 4 byte value
    }
    BOOST_LOCALE_START_CONST_CONDITION
    if(sizeof(long) > 4) {
        BOOST_LOCALE_END_CONST_CONDITION
        TEST_EQ_EXPR("4294967295", 4294967295); // largest 4 byte value
        TEST_EQ_EXPR("100000000000", 100000000000);
    }

    // Unary
    TEST_EQ_EXPR("!0", 1);
    TEST_EQ_EXPR("!1", 0);
    TEST_EQ_EXPR("!42", 0);
    // Arithmetic
    TEST_EQ_EXPR("30 / 5", 6);
    TEST_EQ_EXPR("3 * 7", 21);
    TEST_EQ_EXPR("14 % 3", 2);
    TEST_EQ_EXPR("2 + 5", 7);
    TEST_EQ_EXPR("7 - 5", 2);
    // Comparison
    TEST_EQ_EXPR("30 > 5", 1);
    TEST_EQ_EXPR("5 > 5", 0);
    TEST_EQ_EXPR("4 > 5", 0);
    TEST_EQ_EXPR("30 < 5", 0);
    TEST_EQ_EXPR("5 < 5", 0);
    TEST_EQ_EXPR("4 < 5", 1);
    TEST_EQ_EXPR("30 >= 5", 1);
    TEST_EQ_EXPR("5 >= 5", 1);
    TEST_EQ_EXPR("4 >= 5", 0);
    TEST_EQ_EXPR("30 <= 5", 0);
    TEST_EQ_EXPR("5 <= 5", 1);
    TEST_EQ_EXPR("4 <= 5", 1);
    TEST_EQ_EXPR("5 == 5", 1);
    TEST_EQ_EXPR("4 == 5", 0);
    TEST_EQ_EXPR("5 != 5", 0);
    TEST_EQ_EXPR("4 != 5", 1);
    // Logicals
    TEST_EQ_EXPR("0 || 0", 0);
    TEST_EQ_EXPR("0 || 1", 1);
    TEST_EQ_EXPR("1 || 0", 1);
    TEST_EQ_EXPR("1 || 1", 1);
    TEST_EQ_EXPR("0 && 0", 0);
    TEST_EQ_EXPR("0 && 1", 0);
    TEST_EQ_EXPR("1 && 0", 0);
    TEST_EQ_EXPR("1 && 1", 1);
    // Corner case for div/mod
    TEST_EQ_EXPR("0 / 0", 0);
    TEST_EQ_EXPR("1 / 0", 0);
    TEST_EQ_EXPR("20 / 0", 0);
    TEST_EQ_EXPR("0 % 0", 0);
    TEST_EQ_EXPR("1 % 0", 0);
    TEST_EQ_EXPR("20 % 0", 0);
    // Operator precedence
    {
        // Unary over all
        TEST_EQ_EXPR("5 * -3", -15);
        TEST_EQ_EXPR("2 + -3", -1);
        TEST_EQ_EXPR("-(1-6) % --2", 1);
        TEST_EQ_EXPR("!-1", 0);
        TEST_EQ_EXPR("-!0", -1);
        TEST_EQ_EXPR("-!1", 0);
        TEST_EQ_EXPR("5 * !3", 0);
        TEST_EQ_EXPR("5 * !0", 5);
        TEST_EQ_EXPR("2 + !3", 2);
        TEST_EQ_EXPR("2 + !0", 3);
        TEST_EQ_EXPR("!(1-6) % 10", 0);
        TEST_EQ_EXPR("!!(1-6) % 10", 1);
        TEST_EQ_EXPR("!0 == 1", 1);
        TEST_EQ_EXPR("!0 || 1", 1);
        // Mul over add/sub
        TEST_EQ_EXPR("5 * 3 + 1", 16);
        TEST_EQ_EXPR("1 + 3 * 5", 16);
        TEST_EQ_EXPR("5 * 3 - 1", 14);
        TEST_EQ_EXPR("10 - 2 * 3", 4);
        // Div over add/sub
        TEST_EQ_EXPR("30 / 2 + 1", 16);
        TEST_EQ_EXPR("1 + 30 / 2", 16);
        TEST_EQ_EXPR("30 / 2 - 1", 14);
        TEST_EQ_EXPR("10 - 30 / 2", -5);
        // Mod over add/sub
        TEST_EQ_EXPR("5 % 3 + 1", 3);
        TEST_EQ_EXPR("1 + 5 % 3", 3);
        TEST_EQ_EXPR("5 % 4 - 1", 0);
        TEST_EQ_EXPR("5 - 5 % 3", 3);
        // Same precedence of add/sub
        TEST_EQ_EXPR("10 - 4 - 2", 4);
        // Same precedence of mul/div/mod
        TEST_EQ_EXPR("2 * 20 / 40", 1);
        TEST_EQ_EXPR("20 / 20 * 2", 2);
        TEST_EQ_EXPR("7 * 2 % 10", 4);
        TEST_EQ_EXPR("7 % 4 * 2", 6);
        // Comparison after operation
        TEST_EQ_EXPR("5 * 2 == 5 + 5", 1);
        TEST_EQ_EXPR("50 / 5 != 40 % 40", 1);
        TEST_EQ_EXPR("5 * 2 < 5 + 6", 1);
        TEST_EQ_EXPR("5 * 2 <= 5 + 6", 1);
        TEST_EQ_EXPR("50 / 5 > 40 % 40", 1);
        TEST_EQ_EXPR("50 / 5 >= 40 % 40", 1);
        // Relational before equal with counter example
        TEST_EQ_EXPR("1 < 2 == 1", 1);
        TEST_EQ_EXPR("1 < (2 == 1)", 0);
        TEST_EQ_EXPR("0 == 0 < 0", 1);
        TEST_EQ_EXPR("(0 == 0) < 0", 0);

        TEST_EQ_EXPR("1 <= 2 == 1", 1);
        TEST_EQ_EXPR("1 <= (2 == 1)", 0);
        TEST_EQ_EXPR("0 == 0 <= -1", 1);
        TEST_EQ_EXPR("(0 == 0) <= -1", 0);

        TEST_EQ_EXPR("0 > 0 == 0", 1);
        TEST_EQ_EXPR("0 > (0 == 0)", 0);
        TEST_EQ_EXPR("1 == 2 > 1", 1);
        TEST_EQ_EXPR("(1 == 2) > 1", 0);

        TEST_EQ_EXPR("-1 >= 0 == 0", 1);
        TEST_EQ_EXPR("-1 >= (0 == 0)", 0);
        TEST_EQ_EXPR("1 == 2 >= 1", 1);
        TEST_EQ_EXPR("(1 == 2) >= 1", 0);
        // Comparison on both sides
        TEST_EQ_EXPR("0 > -2 == -1 < 0", 1);
        TEST_EQ_EXPR("((0 > -2) == -1) < 0", 0); // RHS late
        TEST_EQ_EXPR("0 > (-2 == (-1 < 0))", 0); // LHS late
        TEST_EQ_EXPR("(0 > (-2 == -1)) < 0", 0); // both late
        // Logical after comparison with counter example
        TEST_EQ_EXPR("1 && 3==3", 1);
        TEST_EQ_EXPR("(1 && 3)==3", 0);
        TEST_EQ_EXPR("3==3 && 1", 1);
        TEST_EQ_EXPR("3==(3 && 1)", 0);
        TEST_EQ_EXPR("0 || 3==3", 1);
        TEST_EQ_EXPR("(0 || 3)==3", 0);
        TEST_EQ_EXPR("3==3 || 0", 1);
        TEST_EQ_EXPR("3==(3 || 0)", 0);
        // OR after AND
        TEST_EQ_EXPR("0 && 0 || 1", 1);
        TEST_EQ_EXPR("0 && (0 || 1)", 0);
        TEST_EQ_EXPR("1 || 0 && 0", 1);
    }
#undef TEST_EQ_EXPR

    // Random test using the variable comparing against C++ evaluated result
#define TEST_PLURAL_EXPR(expr) \
    test_plural_expr_rand(     \
      [](long long n) {        \
          (void)n;             \
          return expr;         \
      },                       \
      #expr);
    TEST_PLURAL_EXPR(42);
    TEST_PLURAL_EXPR(1337);
    TEST_PLURAL_EXPR(n + 3);
    TEST_PLURAL_EXPR(n - 3);
    TEST_PLURAL_EXPR(n * 5);
    TEST_PLURAL_EXPR(n / 5);
    TEST_PLURAL_EXPR(n % 8);
    // Parentheses
    TEST_PLURAL_EXPR((5 * n) + 3);
    TEST_PLURAL_EXPR(5 * (n + 3));
    // Comparisons and ternary
    TEST_PLURAL_EXPR(n % 2 == 0 ? n + 5 : n * 3);
    TEST_PLURAL_EXPR(n % 2 != 0 ? n + 5 : n * 3);
    TEST_PLURAL_EXPR(n % 4 < 2 ? n + 5 : n * 3);
    TEST_PLURAL_EXPR(n % 4 <= 2 ? n + 5 : n * 3);
    TEST_PLURAL_EXPR(n % 4 > 2 ? n + 5 : n * 3);
    TEST_PLURAL_EXPR(n % 4 >= 2 ? n + 5 : n * 3);
    // Complex expression (e.g. for Russian)
    TEST_PLURAL_EXPR((n % 10 == 1 && n % 100 != 11                                  ? 0 :
                      n % 10 >= 2 && n % 10 <= 4 && (n % 100 < 10 || n % 100 >= 20) ? 1 :
                                                                                      2));
#undef TEST_PLURAL_EXPR

    constexpr auto minVal = std::numeric_limits<long long>::min();
    constexpr auto maxVal = std::numeric_limits<long long>::max();

    // E.g. Japanese
    {
        const auto p = COMPILE_PLURAL_EXPR("0");
        TEST_EQ(p(0), 0);
        TEST_EQ(p(minVal), 0);
        TEST_EQ(p(maxVal), 0);
    }
    // E.g. English
    {
        const auto p = COMPILE_PLURAL_EXPR("(n != 1)");
        TEST_EQ(p(0), 1);
        TEST_EQ(p(1), 0);
        TEST_EQ(p(minVal), 1);
        TEST_EQ(p(maxVal), 1);
    }
    // E.g. French
    {
        const auto p = COMPILE_PLURAL_EXPR("(n > 1)");
        TEST_EQ(p(0), 0);
        TEST_EQ(p(1), 0);
        TEST_EQ(p(2), 1);
        TEST_EQ(p(minVal), 0);
        TEST_EQ(p(maxVal), 1);
    }
    // E.g. Latvian
    {
        const auto p = COMPILE_PLURAL_EXPR("(n%10==1 && n%100!=11 ? 0 : n != 0 ? 1 : 2)");
        TEST_EQ(p(0), 2);
        TEST_EQ(p(1), 0);
        TEST_EQ(p(2), 1);
        TEST_EQ(p(3), 1);
        TEST_EQ(p(11), 1);
        TEST_EQ(p(12), 1);
        TEST_EQ(p(21), 0);
        TEST_EQ(p(31), 0);
        TEST_EQ(p(101), 0);
        TEST_EQ(p(111), 1);
    }
    // E.g. Irish
    {
        const auto p = COMPILE_PLURAL_EXPR("n == 1 ? 0 : n == 2 ? 1 : 2");
        TEST_EQ(p(0), 2);
        TEST_EQ(p(1), 0);
        TEST_EQ(p(2), 1);
        TEST_EQ(p(3), 2);
        TEST_EQ(p(4), 2);
        TEST_EQ(p(minVal), 2);
        TEST_EQ(p(maxVal), 2);
    }
    // E.g. Czech
    {
        const auto p = COMPILE_PLURAL_EXPR("(n==1) ? 0 : (n>=2 && n<=4) ? 1 : 2");
        TEST_EQ(p(0), 2);
        TEST_EQ(p(1), 0);
        TEST_EQ(p(2), 1);
        TEST_EQ(p(3), 1);
        TEST_EQ(p(4), 1);
        TEST_EQ(p(5), 2);
        TEST_EQ(p(minVal), 2);
        TEST_EQ(p(maxVal), 2);
    }
#undef COMPILE_PLURAL_EXPR
    // Error cases
    TEST(!compile("") && compile("n")); // Empty
    // Invalid comparison
    TEST(!compile("n===1") && compile("n==1"));
    TEST(!compile("n!1") && compile("n!=1"));
    TEST(!compile("n!<1") && compile("n<1"));
    TEST(!compile("n<==1") && compile("n<=1"));
    TEST(!compile("n<>1") && compile("n>1"));
    // Incomplete ternary
    TEST(!compile("n==1 ?"));
    TEST(!compile("n==1 ? 1"));
    TEST(!compile("n==1 ? 1 :"));
    TEST(compile("n==1 ? 1 : 0"));
    // Missing closing parenthesis
    TEST(!compile("(n==1") && compile("(n==1)"));
    TEST(!compile("(n + 1") && compile("(n + 1)"));
    TEST(!compile("(n==1 ? 1 : 2") && compile("(n==1 ? 1 : 2)"));
    // Extra closing parenthesis
    TEST(!compile("n==1)") && compile("(n==1)"));
    TEST(!compile("n + 1)") && compile("(n + 1)"));
    TEST(!compile("n==1 ? 1 : 2)") && compile("(n==1 ? 1 : 2)"));
    // Empty parenthesis
    TEST(!compile("n==()1"));
    TEST(!compile("n==()"));
    // Missing operator for unary op
    TEST(!compile("1==!"));
    TEST(!compile("1 + -"));
    // No bitwise ops
    TEST(!compile("n << 1"));
    TEST(!compile("n >> 1"));
    TEST(!compile("n & 1"));
    TEST(!compile("n | 1"));
    TEST(!compile("n ^ 1"));
    TEST(!compile("~n"));
}

void test_main(int /*argc*/, char** /*argv*/)
{
    test_plural_expr();
}
