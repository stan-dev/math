#include <stan/math/rev/scal.hpp>
#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <string>

////////////////////////////////
/**
 * Require var
 */
template <typename T1, typename = void>
struct require_var_tester : std::false_type {};

template <typename T1>
struct require_var_tester<T1, stan::require_var<T1>> : std::true_type {};

TEST(requires, var_test) {
  using stan::math::var;
  EXPECT_FALSE((require_var_tester<std::string>::value));
  EXPECT_FALSE((require_var_tester<double>::value));
  EXPECT_FALSE((require_var_tester<int>::value));
  EXPECT_TRUE((require_var_tester<var>::value));
}

/**
 * Require not var
 */
template <typename T1, typename = void>
struct require_not_var_tester : std::false_type {};

template <typename T1>
struct require_not_var_tester<T1, stan::require_not_var<T1>> : std::true_type {
};

TEST(requires, not_var_test) {
  using stan::math::var;
  EXPECT_TRUE((require_not_var_tester<std::string>::value));
  EXPECT_TRUE((require_not_var_tester<double>::value));
  EXPECT_TRUE((require_not_var_tester<int>::value));
  EXPECT_FALSE((require_not_var_tester<var>::value));
}

/**
 * Require all var
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_all_var_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_all_var_tester<T1, T2, T3, stan::require_all_var<T1, T2, T3>>
    : std::true_type {};

TEST(requires, all_var_test) {
  using stan::math::var;
  EXPECT_FALSE((require_all_var_tester<double, std::string, double>::value));
  EXPECT_FALSE((require_all_var_tester<double, double, double>::value));
  EXPECT_FALSE((require_all_var_tester<int, int, int>::value));
  EXPECT_TRUE((require_all_var_tester<var, var, var>::value));
}

/**
 * Require Not all var
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_all_not_var_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_all_not_var_tester<T1, T2, T3,
                                  stan::require_all_not_var<T1, T2, T3>>
    : std::true_type {};

TEST(requires, not_all_var_test) {
  using stan::math::var;
  EXPECT_TRUE((require_all_not_var_tester<var, std::string, double>::value));
  EXPECT_TRUE((require_all_not_var_tester<double, var, double>::value));
  EXPECT_TRUE((require_all_not_var_tester<int, int, var>::value));
  EXPECT_FALSE((require_all_not_var_tester<var, var, var>::value));
}

/**
 * Require any var
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_any_var_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_any_var_tester<T1, T2, T3, stan::require_any_var<T1, T2, T3>>
    : std::true_type {};

TEST(requires, any_var_test) {
  using stan::math::var;
  EXPECT_FALSE(
      (require_any_var_tester<std::string, std::string, std::string>::value));
  EXPECT_FALSE((require_any_var_tester<double, std::string, double>::value));
  EXPECT_FALSE((require_any_var_tester<int, int, int>::value));
  EXPECT_TRUE((require_any_var_tester<var, int, int>::value));
  EXPECT_TRUE((require_any_var_tester<var, int, var>::value));
}

/**
 * Require Not any var
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_any_not_var_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_any_not_var_tester<T1, T2, T3,
                                  stan::require_any_not_var<T1, T2, T3>>
    : std::true_type {};

TEST(requires, not_any_var_test) {
  using stan::math::var;
  EXPECT_TRUE((require_any_not_var_tester<std::string, std::string,
                                          std::string>::value));
  EXPECT_TRUE((require_any_not_var_tester<double, double, double>::value));
  EXPECT_TRUE((require_any_not_var_tester<int, int, int>::value));
  EXPECT_FALSE((require_any_not_var_tester<int, var, var>::value));
}

////////////////////////////////
/**
 * Require var
 */
template <typename T1, typename = void>
struct require_var_or_fvar_tester : std::false_type {};

template <typename T1>
struct require_var_or_fvar_tester<T1, stan::require_var_or_fvar<T1>>
    : std::true_type {};

TEST(requires, var_or_fvar_test) {
  using stan::math::var;
  EXPECT_FALSE((require_var_or_fvar_tester<std::string>::value));
  EXPECT_FALSE((require_var_or_fvar_tester<double>::value));
  EXPECT_FALSE((require_var_or_fvar_tester<int>::value));
  EXPECT_TRUE((require_var_or_fvar_tester<var>::value));
}

/**
 * Require not var
 */
template <typename T1, typename = void>
struct require_not_var_or_fvar_tester : std::false_type {};

template <typename T1>
struct require_not_var_or_fvar_tester<T1, stan::require_not_var_or_fvar<T1>>
    : std::true_type {};

TEST(requires, not_var_or_fvar_test) {
  using stan::math::var;
  EXPECT_TRUE((require_not_var_or_fvar_tester<std::string>::value));
  EXPECT_TRUE((require_not_var_or_fvar_tester<double>::value));
  EXPECT_TRUE((require_not_var_or_fvar_tester<int>::value));
  EXPECT_FALSE((require_not_var_or_fvar_tester<var>::value));
}

/**
 * Require all var
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_all_var_or_fvar_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_all_var_or_fvar_tester<T1, T2, T3,
                                      stan::require_all_var_or_fvar<T1, T2, T3>>
    : std::true_type {};

TEST(requires, all_var_or_fvar_test) {
  using stan::math::var;
  EXPECT_FALSE(
      (require_all_var_or_fvar_tester<double, std::string, double>::value));
  EXPECT_FALSE((require_all_var_or_fvar_tester<double, double, double>::value));
  EXPECT_FALSE((require_all_var_or_fvar_tester<int, int, int>::value));
  EXPECT_TRUE((require_all_var_or_fvar_tester<var, var, var>::value));
}

/**
 * Require Not all var
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_all_not_var_or_fvar_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_all_not_var_or_fvar_tester<
    T1, T2, T3, stan::require_all_not_var_or_fvar<T1, T2, T3>>
    : std::true_type {};

TEST(requires, not_all_var_or_fvar_test) {
  using stan::math::var;
  EXPECT_TRUE(
      (require_all_not_var_or_fvar_tester<var, std::string, double>::value));
  EXPECT_TRUE((require_all_not_var_or_fvar_tester<double, var, double>::value));
  EXPECT_TRUE((require_all_not_var_or_fvar_tester<int, int, var>::value));
  EXPECT_FALSE((require_all_not_var_or_fvar_tester<var, var, var>::value));
}

/**
 * Require any var
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_any_var_or_fvar_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_any_var_or_fvar_tester<T1, T2, T3,
                                      stan::require_any_var_or_fvar<T1, T2, T3>>
    : std::true_type {};

TEST(requires, any_var_or_fvar_test) {
  using stan::math::var;
  EXPECT_FALSE((require_any_var_or_fvar_tester<std::string, std::string,
                                               std::string>::value));
  EXPECT_FALSE(
      (require_any_var_or_fvar_tester<double, std::string, double>::value));
  EXPECT_FALSE((require_any_var_or_fvar_tester<int, int, int>::value));
  EXPECT_TRUE((require_any_var_or_fvar_tester<var, int, int>::value));
  EXPECT_TRUE((require_any_var_or_fvar_tester<var, int, var>::value));
}

/**
 * Require Not any var
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_any_not_var_or_fvar_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_any_not_var_or_fvar_tester<
    T1, T2, T3, stan::require_any_not_var_or_fvar<T1, T2, T3>>
    : std::true_type {};

TEST(requires, not_any_var_or_fvar_test) {
  using stan::math::var;
  EXPECT_TRUE((require_any_not_var_or_fvar_tester<std::string, std::string,
                                                  std::string>::value));
  EXPECT_TRUE(
      (require_any_not_var_or_fvar_tester<double, double, double>::value));
  EXPECT_TRUE((require_any_not_var_or_fvar_tester<int, int, int>::value));
  EXPECT_FALSE((require_any_not_var_or_fvar_tester<int, var, var>::value));
}

////////////////////////////////
/**
 * Require var
 */
template <typename T1, typename = void>
struct require_stan_scalar_tester : std::false_type {};

template <typename T1>
struct require_stan_scalar_tester<T1, stan::require_stan_scalar<T1>>
    : std::true_type {};

TEST(requires, stan_scalar_test) {
  using stan::math::var;
  EXPECT_FALSE((require_stan_scalar_tester<std::string>::value));
  EXPECT_TRUE((require_stan_scalar_tester<double>::value));
  EXPECT_TRUE((require_stan_scalar_tester<int>::value));
  EXPECT_TRUE((require_stan_scalar_tester<var>::value));
}

/**
 * Require not var
 */
template <typename T1, typename = void>
struct require_not_stan_scalar_tester : std::false_type {};

template <typename T1>
struct require_not_stan_scalar_tester<T1, stan::require_not_stan_scalar<T1>>
    : std::true_type {};

TEST(requires, not_stan_scalar_test) {
  using stan::math::var;
  EXPECT_TRUE((require_not_stan_scalar_tester<std::string>::value));
  EXPECT_FALSE((require_not_stan_scalar_tester<double>::value));
  EXPECT_FALSE((require_not_stan_scalar_tester<int>::value));
  EXPECT_FALSE((require_not_stan_scalar_tester<var>::value));
}

/**
 * Require all var
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_all_stan_scalar_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_all_stan_scalar_tester<T1, T2, T3,
                                      stan::require_all_stan_scalar<T1, T2, T3>>
    : std::true_type {};

TEST(requires, all_stan_scalar_test) {
  using stan::math::var;
  EXPECT_FALSE(
      (require_all_stan_scalar_tester<double, std::string, double>::value));
  EXPECT_TRUE((require_all_stan_scalar_tester<double, double, double>::value));
  EXPECT_TRUE((require_all_stan_scalar_tester<int, int, int>::value));
  EXPECT_TRUE((require_all_stan_scalar_tester<var, var, var>::value));
}

/**
 * Require Not all var
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_all_not_stan_scalar_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_all_not_stan_scalar_tester<
    T1, T2, T3, stan::require_all_not_stan_scalar<T1, T2, T3>>
    : std::true_type {};

TEST(requires, not_all_stan_scalar_test) {
  using stan::math::var;
  EXPECT_TRUE(
      (require_all_not_stan_scalar_tester<var, std::string, double>::value));
  EXPECT_FALSE(
      (require_all_not_stan_scalar_tester<double, var, double>::value));
  EXPECT_FALSE((require_all_not_stan_scalar_tester<int, int, var>::value));
  EXPECT_FALSE((require_all_not_stan_scalar_tester<var, var, var>::value));
}

/**
 * Require any var
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_any_stan_scalar_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_any_stan_scalar_tester<T1, T2, T3,
                                      stan::require_any_stan_scalar<T1, T2, T3>>
    : std::true_type {};

TEST(requires, any_stan_scalar_test) {
  using stan::math::var;
  EXPECT_FALSE((require_any_stan_scalar_tester<std::string, std::string,
                                               std::string>::value));
  EXPECT_TRUE(
      (require_any_stan_scalar_tester<double, std::string, double>::value));
  EXPECT_TRUE((require_any_stan_scalar_tester<int, int, int>::value));
  EXPECT_TRUE((require_any_stan_scalar_tester<var, int, int>::value));
  EXPECT_TRUE((require_any_stan_scalar_tester<var, int, var>::value));
}

/**
 * Require Not any var
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_any_not_stan_scalar_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_any_not_stan_scalar_tester<
    T1, T2, T3, stan::require_any_not_stan_scalar<T1, T2, T3>>
    : std::true_type {};

TEST(requires, not_any_stan_scalar_test) {
  using stan::math::var;
  EXPECT_TRUE((require_any_not_stan_scalar_tester<std::string, std::string,
                                                  std::string>::value));
  EXPECT_FALSE(
      (require_any_not_stan_scalar_tester<double, double, double>::value));
  EXPECT_FALSE((require_any_not_stan_scalar_tester<int, int, int>::value));
  EXPECT_FALSE((require_any_not_stan_scalar_tester<int, var, var>::value));
}
