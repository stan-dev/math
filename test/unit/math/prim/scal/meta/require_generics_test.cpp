#include <stan/math/prim/scal/meta/require_generics.hpp>
#include <test/unit/math/require_util.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <string>

// Basic tests for the underlying requires
TEST(requires, t_test) {
  using stan::require_t;
  using stan::test::unary_require_tester;
  EXPECT_TRUE((unary_require_tester<require_t, std::true_type>::value));
  EXPECT_FALSE((unary_require_tester<require_t, std::false_type>::value));
}

TEST(requires, not_test) {
  using stan::require_not_t;
  using stan::test::unary_require_tester;
  EXPECT_TRUE((unary_require_tester<require_not_t, std::false_type>::value));
  EXPECT_FALSE((unary_require_tester<require_not_t, std::true_type>::value));
}

TEST(requires, all_not_test) {
  using stan::require_all_not_t;
  using stan::test::variadic_require_tester;
  EXPECT_TRUE((variadic_require_tester<require_all_not_t, std::false_type,
                                       std::false_type>::value));
  EXPECT_FALSE((variadic_require_tester<require_all_not_t, std::true_type,
                                        std::true_type>::value));
  EXPECT_TRUE((variadic_require_tester<require_all_not_t, std::true_type,
                                       std::false_type>::value));
}

TEST(requires, any_not_test) {
  using stan::require_any_not_t;
  using stan::test::variadic_require_tester;
  EXPECT_TRUE((variadic_require_tester<require_any_not_t, std::false_type,
                                       std::false_type>::value));
  EXPECT_FALSE((variadic_require_tester<require_any_not_t, std::true_type,
                                        std::true_type>::value));
  EXPECT_FALSE((variadic_require_tester<require_any_not_t, std::true_type,
                                        std::false_type>::value));
}

TEST(requires, all_test) {
  using stan::require_all_t;
  using stan::test::variadic_require_tester;
  EXPECT_FALSE((variadic_require_tester<require_all_t, std::false_type,
                                        std::false_type>::value));
  EXPECT_TRUE((variadic_require_tester<require_all_t, std::true_type,
                                       std::true_type>::value));
  EXPECT_FALSE((variadic_require_tester<require_all_t, std::true_type,
                                        std::false_type>::value));
}

TEST(requires, any_test) {
  using stan::require_any_t;
  using stan::test::variadic_require_tester;
  EXPECT_FALSE((variadic_require_tester<require_any_t, std::false_type,
                                        std::false_type>::value));
  EXPECT_TRUE((variadic_require_tester<require_any_t, std::true_type,
                                       std::true_type>::value));
  EXPECT_TRUE((variadic_require_tester<require_any_t, std::true_type,
                                       std::false_type>::value));
}

// Test same
TEST(requires, same_test) {
  using stan::require_same_t;
  using stan::test::variadic_require_tester;
  EXPECT_FALSE(
      (variadic_require_tester<stan::require_same_t, double, int>::value));
  EXPECT_TRUE(
      (variadic_require_tester<stan::require_same_t, double, double>::value));
  EXPECT_TRUE((variadic_require_tester<stan::require_same_t, int, int>::value));
}

TEST(requires, not_same_test) {
  using stan::require_not_same_t;
  using stan::test::variadic_require_tester;
  EXPECT_TRUE(
      (variadic_require_tester<require_not_same_t, double, int>::value));
  EXPECT_FALSE(
      (variadic_require_tester<require_not_same_t, double, double>::value));
  EXPECT_FALSE((variadic_require_tester<require_not_same_t, int, int>::value));
}

TEST(requires, all_same_test) {
  using stan::require_all_same_t;
  using stan::test::variadic_require_tester;
  EXPECT_FALSE((variadic_require_tester<require_all_same_t, double, std::string,
                                        double>::value));
  EXPECT_TRUE((variadic_require_tester<require_all_same_t, double, double,
                                       double>::value));
  EXPECT_TRUE(
      (variadic_require_tester<require_all_same_t, int, int, int>::value));
}

TEST(requires, all_not_same_test) {
  using stan::require_all_not_same_t;
  using stan::test::variadic_require_tester;
  EXPECT_TRUE((variadic_require_tester<require_all_not_same_t, double, int,
                                       double>::value));
  EXPECT_FALSE((variadic_require_tester<require_all_not_same_t, double, double,
                                        double>::value));
  EXPECT_FALSE(
      (variadic_require_tester<require_all_not_same_t, int, int, int>::value));
}

// Double or Int
TEST(requires, double_or_int_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_double_or_int_t>::unary();
  require_scal_checker<stan::require_not_double_or_int_t>::not_unary();
  require_scal_checker<stan::require_all_double_or_int_t>::all();
  require_scal_checker<stan::require_all_not_double_or_int_t>::all_not();
  require_scal_checker<stan::require_any_double_or_int_t>::any();
  require_scal_checker<stan::require_any_not_double_or_int_t>::any_not();
}

TEST(requires, arithmetic_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_arithmetic_t>::unary();
  require_scal_checker<stan::require_not_arithmetic_t>::not_unary();
  require_scal_checker<stan::require_all_arithmetic_t>::all();
  require_scal_checker<stan::require_all_not_arithmetic_t>::all_not();
  require_scal_checker<stan::require_any_arithmetic_t>::any();
  require_scal_checker<stan::require_any_not_arithmetic_t>::any_not();
}

TEST(requires, var_or_arithmetic_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_var_or_arithmetic_t>::unary();
  require_scal_checker<stan::require_not_var_or_arithmetic_t>::not_unary();
  require_scal_checker<stan::require_all_var_or_arithmetic_t>::all();
  require_scal_checker<stan::require_all_not_var_or_arithmetic_t>::all_not();
  require_scal_checker<stan::require_any_var_or_arithmetic_t>::any();
  require_scal_checker<stan::require_any_not_var_or_arithmetic_t>::any_not();
}
