#include <stan/math/fwd/scal.hpp>
#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <string>

////////////////////////////////
/**
 * Require fvar
 */
template <typename T1, typename = void>
struct require_fvar_tester : std::false_type {};

template <typename T1>
struct require_fvar_tester<T1, stan::require_fvar<T1>> : std::true_type {};

TEST(requires, fvar_test) {
  using stan::math::fvar;
  EXPECT_FALSE((require_fvar_tester<std::string>::value));
  EXPECT_FALSE((require_fvar_tester<double>::value));
  EXPECT_FALSE((require_fvar_tester<int>::value));
  EXPECT_TRUE((require_fvar_tester<fvar<double>>::value));
}

/**
 * Require not fvar
 */
template <typename T1, typename = void>
struct require_not_fvar_tester : std::false_type {};

template <typename T1>
struct require_not_fvar_tester<T1, stan::require_not_fvar<T1>>
    : std::true_type {};

TEST(requires, not_fvar_test) {
  using stan::math::fvar;
  EXPECT_TRUE((require_not_fvar_tester<std::string>::value));
  EXPECT_TRUE((require_not_fvar_tester<double>::value));
  EXPECT_TRUE((require_not_fvar_tester<int>::value));
  EXPECT_FALSE((require_not_fvar_tester<fvar<double>>::value));
}

/**
 * Require all fvar
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_all_fvar_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_all_fvar_tester<T1, T2, T3, stan::require_all_fvar<T1, T2, T3>>
    : std::true_type {};

TEST(requires, all_fvar_test) {
  using stan::math::fvar;
  EXPECT_FALSE((require_all_fvar_tester<double, std::string, double>::value));
  EXPECT_FALSE((require_all_fvar_tester<double, double, double>::value));
  EXPECT_FALSE((require_all_fvar_tester<int, int, int>::value));
  EXPECT_TRUE((require_all_fvar_tester<fvar<double>, fvar<double>,
                                       fvar<double>>::value));
}

/**
 * Require Not all fvar
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_all_not_fvar_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_all_not_fvar_tester<T1, T2, T3,
                                   stan::require_all_not_fvar<T1, T2, T3>>
    : std::true_type {};

TEST(requires, not_all_fvar_test) {
  using stan::math::fvar;
  EXPECT_TRUE(
      (require_all_not_fvar_tester<fvar<double>, std::string, double>::value));
  EXPECT_TRUE(
      (require_all_not_fvar_tester<double, fvar<double>, double>::value));
  EXPECT_TRUE((require_all_not_fvar_tester<int, int, fvar<double>>::value));
  EXPECT_FALSE((require_all_not_fvar_tester<fvar<double>, fvar<double>,
                                            fvar<double>>::value));
}

/**
 * Require any fvar
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_any_fvar_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_any_fvar_tester<T1, T2, T3, stan::require_any_fvar<T1, T2, T3>>
    : std::true_type {};

TEST(requires, any_fvar_test) {
  using stan::math::fvar;
  EXPECT_FALSE(
      (require_any_fvar_tester<std::string, std::string, std::string>::value));
  EXPECT_FALSE((require_any_fvar_tester<double, std::string, double>::value));
  EXPECT_FALSE((require_any_fvar_tester<int, int, int>::value));
  EXPECT_TRUE((require_any_fvar_tester<fvar<double>, int, int>::value));
  EXPECT_TRUE(
      (require_any_fvar_tester<fvar<double>, int, fvar<double>>::value));
}

/**
 * Require Not any fvar
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_any_not_fvar_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_any_not_fvar_tester<T1, T2, T3,
                                   stan::require_any_not_fvar<T1, T2, T3>>
    : std::true_type {};

TEST(requires, not_any_fvar_test) {
  using stan::math::fvar;
  EXPECT_TRUE((require_any_not_fvar_tester<std::string, std::string,
                                           std::string>::value));
  EXPECT_TRUE((require_any_not_fvar_tester<double, double, double>::value));
  EXPECT_TRUE((require_any_not_fvar_tester<int, int, int>::value));
  EXPECT_FALSE(
      (require_any_not_fvar_tester<int, fvar<double>, fvar<double>>::value));
}
