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
  using stan::test::require_variadic_checker;
  EXPECT_TRUE((require_variadic_checker<require_all_not_t, std::false_type,
                                        std::false_type>::value));
  EXPECT_FALSE((require_variadic_checker<require_all_not_t, std::true_type,
                                         std::true_type>::value));
  EXPECT_TRUE((require_variadic_checker<require_all_not_t, std::true_type,
                                        std::false_type>::value));
}

TEST(requires, any_not_test) {
  using stan::require_any_not_t;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE((require_variadic_checker<require_any_not_t, std::false_type,
                                        std::false_type>::value));
  EXPECT_FALSE((require_variadic_checker<require_any_not_t, std::true_type,
                                         std::true_type>::value));
  EXPECT_FALSE((require_variadic_checker<require_any_not_t, std::true_type,
                                         std::false_type>::value));
}

TEST(requires, all_test) {
  using stan::require_all_t;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE((require_variadic_checker<require_all_t, std::false_type,
                                         std::false_type>::value));
  EXPECT_TRUE((require_variadic_checker<require_all_t, std::true_type,
                                        std::true_type>::value));
  EXPECT_FALSE((require_variadic_checker<require_all_t, std::true_type,
                                         std::false_type>::value));
}

TEST(requires, any_test) {
  using stan::require_any_t;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE((require_variadic_checker<require_any_t, std::false_type,
                                         std::false_type>::value));
  EXPECT_TRUE((require_variadic_checker<require_any_t, std::true_type,
                                        std::true_type>::value));
  EXPECT_TRUE((require_variadic_checker<require_any_t, std::true_type,
                                        std::false_type>::value));
}

// Test same
TEST(requires, same_test) {
  using stan::require_same_t;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE(
      (require_variadic_checker<stan::require_same_t, double, int>::value));
  EXPECT_TRUE(
      (require_variadic_checker<stan::require_same_t, double, double>::value));
  EXPECT_TRUE(
      (require_variadic_checker<stan::require_same_t, int, int>::value));
}

TEST(requires, not_same_test) {
  using stan::require_not_same_t;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE(
      (require_variadic_checker<require_not_same_t, double, int>::value));
  EXPECT_FALSE(
      (require_variadic_checker<require_not_same_t, double, double>::value));
  EXPECT_FALSE((require_variadic_checker<require_not_same_t, int, int>::value));
}

TEST(requires, all_same_test) {
  using stan::require_all_same_t;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE((require_variadic_checker<require_all_same_t, double,
                                         std::string, double>::value));
  EXPECT_TRUE((require_variadic_checker<require_all_same_t, double, double,
                                        double>::value));
  EXPECT_TRUE(
      (require_variadic_checker<require_all_same_t, int, int, int>::value));
}

TEST(requires, all_not_same_test) {
  using stan::require_any_not_same_t;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE((require_variadic_checker<require_any_not_same_t, double, int,
                                        double>::value));
  EXPECT_FALSE((require_variadic_checker<require_any_not_same_t, double, double,
                                         double>::value));
  EXPECT_FALSE(
      (require_variadic_checker<require_any_not_same_t, int, int, int>::value));
}

// Test convertible
TEST(requires, convertible_test) {
  using stan::require_convertible_t;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE((require_variadic_checker<stan::require_convertible_t, double,
                                         char[1]>::value));
  EXPECT_TRUE((require_variadic_checker<stan::require_convertible_t, double,
                                        double>::value));
  EXPECT_TRUE(
      (require_variadic_checker<stan::require_convertible_t, int, int>::value));
}

TEST(requires, not_convertible_test) {
  using stan::require_not_convertible_t;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE((require_variadic_checker<require_not_convertible_t, double,
                                        char[1]>::value));
  EXPECT_FALSE((require_variadic_checker<require_not_convertible_t, double,
                                         double>::value));
  EXPECT_FALSE(
      (require_variadic_checker<require_not_convertible_t, int, int>::value));
}

TEST(requires, all_convertible_test) {
  using stan::require_all_convertible_t;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE((require_variadic_checker<require_all_convertible_t, double,
                                         std::string, double>::value));
  EXPECT_TRUE((require_variadic_checker<require_all_convertible_t, double, int,
                                        int>::value));
  EXPECT_TRUE((require_variadic_checker<require_all_convertible_t, int, int,
                                        int>::value));
}

TEST(requires, all_not_convertible_test) {
  using stan::require_any_not_convertible_t;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE((require_variadic_checker<require_any_not_convertible_t, double,
                                        int, std::string>::value));
  EXPECT_FALSE((require_variadic_checker<require_any_not_convertible_t, double,
                                         double, double>::value));
  EXPECT_FALSE((require_variadic_checker<require_any_not_convertible_t, int,
                                         double, int>::value));
}

// Double or Int
TEST(requires, double_or_int_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_double_or_int_t, double, int>::unary();
  require_scal_checker<stan::require_not_double_or_int_t, double,
                       int>::not_unary();
  require_scal_checker<stan::require_all_double_or_int_t, double, int>::all();
  require_scal_checker<stan::require_all_not_double_or_int_t, double,
                       int>::all_not();
  require_scal_checker<stan::require_any_double_or_int_t, double, int>::any();
  require_scal_checker<stan::require_any_not_double_or_int_t, double,
                       int>::any_not();
}

TEST(requires, arithmetic_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_arithmetic_t, double, float, int>::unary();
  require_scal_checker<stan::require_not_arithmetic_t, double, float,
                       int>::not_unary();
  require_scal_checker<stan::require_all_arithmetic_t, double, float,
                       int>::all();
  require_scal_checker<stan::require_all_not_arithmetic_t, double, float,
                       int>::all_not();
  require_scal_checker<stan::require_any_arithmetic_t, double, float,
                       int>::any();
  require_scal_checker<stan::require_any_not_arithmetic_t, double, float,
                       int>::any_not();
}

TEST(requires, var_or_arithmetic_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_var_or_arithmetic_t, double, int>::unary();
  require_scal_checker<stan::require_not_var_or_arithmetic_t, double,
                       int>::not_unary();
  require_scal_checker<stan::require_all_var_or_arithmetic_t, double,
                       int>::all();
  require_scal_checker<stan::require_all_not_var_or_arithmetic_t, double,
                       int>::all_not();
  require_scal_checker<stan::require_any_var_or_arithmetic_t, double,
                       int>::any();
  require_scal_checker<stan::require_any_not_var_or_arithmetic_t, double,
                       int>::any_not();
}
