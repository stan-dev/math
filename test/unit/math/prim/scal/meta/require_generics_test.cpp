#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <string>

/**
* Require base generic tester
*/
template <template <class> class Check, typename T1, typename = void>
struct unary_require_tester_impl : std::false_type {};

template <template <class> class Check, typename T1>
struct unary_require_tester_impl<Check, T1, decltype((void)(Check<T1>()))> : std::true_type {};

template <template <class> class Check, typename T = void, typename... Types>
struct unary_require_tester {
 static const bool value = (unary_require_tester_impl<Check, T>::value && unary_require_tester<Check, Types...>::value);
};

template <template <class> class Check>
struct unary_require_tester<Check, void> : std::true_type {} ;

TEST(requires, t_test) {
  using stan::require_t;
  EXPECT_TRUE((unary_require_tester<require_t, std::true_type>::value));
  EXPECT_FALSE((unary_require_tester<require_t, std::false_type>::value));
}


TEST(requires, not_test) {
  using stan::require_not_t;
  EXPECT_TRUE((unary_require_tester<require_not_t, std::false_type>::value));
  EXPECT_FALSE((unary_require_tester<require_not_t, std::true_type>::value));
}

/**
 * Require all not
 */
template <typename T1, typename T2, typename = void>
struct require_all_not_tester : std::false_type {};

template <typename T1, typename T2>
struct require_all_not_tester<T1, T2, stan::require_all_not_t<T1, T2>>
    : std::true_type {};

TEST(requires, not_all_test) {
  EXPECT_TRUE(
      (require_all_not_tester<std::false_type, std::false_type>::value));
  EXPECT_FALSE((require_all_not_tester<std::true_type, std::true_type>::value));
  EXPECT_TRUE((require_all_not_tester<std::true_type, std::false_type>::value));
}

/**
 * Require any not
 */
template <typename T1, typename T2, typename = void>
struct require_any_not_tester : std::false_type {};

template <typename T1, typename T2>
struct require_any_not_tester<T1, T2, stan::require_any_not_t<T1, T2>>
    : std::true_type {};

TEST(requires, not_any_test) {
  EXPECT_TRUE(
      (require_any_not_tester<std::false_type, std::false_type>::value));
  EXPECT_FALSE((require_any_not_tester<std::true_type, std::true_type>::value));
  EXPECT_FALSE(
      (require_any_not_tester<std::true_type, std::false_type>::value));
}

/**
 * Require all
 */
template <typename T1, typename T2, typename = void>
struct require_all_tester : std::false_type {};

template <typename T1, typename T2>
struct require_all_tester<T1, T2, stan::require_all_t<T1, T2>>
    : std::true_type {};

TEST(requires, all_test) {
  EXPECT_FALSE((require_all_tester<std::false_type, std::false_type>::value));
  EXPECT_TRUE((require_all_tester<std::true_type, std::true_type>::value));
  EXPECT_FALSE((require_all_tester<std::true_type, std::false_type>::value));
}

/**
 * Require any
 */
template <typename T1, typename T2, typename = void>
struct require_any_tester : std::false_type {};

template <typename T1, typename T2>
struct require_any_tester<T1, T2, stan::require_any_t<T1, T2>>
    : std::true_type {};

TEST(requires, any_test) {
  EXPECT_FALSE((require_any_tester<std::false_type, std::false_type>::value));
  EXPECT_TRUE((require_any_tester<std::true_type, std::true_type>::value));
  EXPECT_TRUE((require_any_tester<std::true_type, std::false_type>::value));
}

/**
 * Require same
 */
template <typename T1, typename T2, typename = void>
struct require_same_tester : std::false_type {};

template <typename T1, typename T2>
struct require_same_tester<T1, T2, stan::require_same<T1, T2>>
    : std::true_type {};


TEST(requires, same_test) {
  EXPECT_FALSE((require_same_tester<double, int>::value));
  EXPECT_TRUE((require_same_tester<double, double>::value));
  EXPECT_TRUE((require_same_tester<int, int>::value));
}

/**
 * Require not same
 */
template <typename T1, typename T2, typename = void>
struct require_not_same_tester : std::false_type {};

template <typename T1, typename T2>
struct require_not_same_tester<T1, T2, stan::require_not_same<T1, T2>>
    : std::true_type {};

TEST(requires, not_same_test) {
  EXPECT_TRUE((require_not_same_tester<double, int>::value));
  EXPECT_FALSE((require_not_same_tester<double, double>::value));
  EXPECT_FALSE((require_not_same_tester<int, int>::value));
}

/**
 * Require all same
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_all_same_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_all_same_tester<T1, T2, T3, stan::require_all_same<T1, T2, T3>>
    : std::true_type {};

TEST(requires, all_same_test) {
  EXPECT_FALSE((require_all_same_tester<double, int, double>::value));
  EXPECT_TRUE((require_all_same_tester<double, double, double>::value));
  EXPECT_TRUE((require_all_same_tester<int, int, int>::value));
}

/**
 * Require Not all same
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_all_not_same_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_all_not_same_tester<T1, T2, T3,
                                   stan::require_all_not_same<T1, T2, T3>>
    : std::true_type {};

TEST(requires, not_all_same_test) {
  EXPECT_TRUE((require_all_not_same_tester<double, int, double>::value));
  EXPECT_FALSE((require_all_not_same_tester<double, double, double>::value));
  EXPECT_FALSE((require_all_not_same_tester<int, int, int>::value));
}

/**
 * Require double_or_int
 */
TEST(requires, double_or_int_test) {
  using stan::require_double_or_int;
  EXPECT_FALSE((unary_require_tester<require_double_or_int, std::string>::value));
  EXPECT_TRUE((unary_require_tester<require_double_or_int, double, int>::value));
}

TEST(requires, not_double_or_int_test) {
  using stan::require_not_double_or_int;
  EXPECT_TRUE((unary_require_tester<require_not_double_or_int, std::string>::value));
  EXPECT_FALSE((unary_require_tester<require_not_double_or_int, double, int>::value));
}

/**
 * Require all double_or_int
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_all_double_or_int_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_all_double_or_int_tester<
    T1, T2, T3, stan::require_all_double_or_int<T1, T2, T3>> : std::true_type {
};

TEST(requires, all_double_or_int_test) {
  EXPECT_FALSE(
      (require_all_double_or_int_tester<double, std::string, double>::value));
  EXPECT_TRUE(
      (require_all_double_or_int_tester<double, double, double>::value));
  EXPECT_TRUE((require_all_double_or_int_tester<int, int, int>::value));
}

/**
 * Require Not all double_or_int
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_all_not_double_or_int_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_all_not_double_or_int_tester<
    T1, T2, T3, stan::require_all_not_double_or_int<T1, T2, T3>>
    : std::true_type {};

TEST(requires, not_all_double_or_int_test) {
  EXPECT_TRUE((require_all_not_double_or_int_tester<double, std::string,
                                                    double>::value));
  EXPECT_FALSE(
      (require_all_not_double_or_int_tester<double, double, double>::value));
  EXPECT_FALSE((require_all_not_double_or_int_tester<int, int, int>::value));
}

/**
 * Require any double_or_int
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_any_double_or_int_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_any_double_or_int_tester<
    T1, T2, T3, stan::require_any_double_or_int<T1, T2, T3>> : std::true_type {
};

TEST(requires, any_double_or_int_test) {
  EXPECT_FALSE((require_any_double_or_int_tester<std::string, std::string,
                                                 std::string>::value));
  EXPECT_TRUE(
      (require_any_double_or_int_tester<double, std::string, double>::value));
  EXPECT_TRUE((require_any_double_or_int_tester<int, int, int>::value));
}

/**
 * Require Not any double_or_int
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_any_not_double_or_int_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_any_not_double_or_int_tester<
    T1, T2, T3, stan::require_any_not_double_or_int<T1, T2, T3>>
    : std::true_type {};

TEST(requires, not_any_double_or_int_test) {
  EXPECT_TRUE((require_any_not_double_or_int_tester<std::string, std::string,
                                                    std::string>::value));
  EXPECT_FALSE(
      (require_any_not_double_or_int_tester<double, double, double>::value));
  EXPECT_FALSE((require_any_not_double_or_int_tester<int, int, int>::value));
}

// Require arithmetic

TEST(requires, arithmetic_test) {
  using stan::require_arithmetic;
  EXPECT_FALSE((unary_require_tester<require_arithmetic, std::string>::value));
  EXPECT_TRUE((unary_require_tester<require_arithmetic, double, int>::value));
}

// Require not arithmetic

TEST(requires, not_arithmetic_test) {
  using stan::require_not_arithmetic;
  EXPECT_TRUE((unary_require_tester<require_not_arithmetic, std::string>::value));
  EXPECT_FALSE((unary_require_tester<require_not_arithmetic, double, int>::value));
}

/**
 * Require all arithmetic
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_all_arithmetic_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_all_arithmetic_tester<T1, T2, T3,
                                     stan::require_all_arithmetic<T1, T2, T3>>
    : std::true_type {};

TEST(requires, all_arithmetic_test) {
  EXPECT_FALSE(
      (require_all_arithmetic_tester<double, std::string, double>::value));
  EXPECT_TRUE((require_all_arithmetic_tester<double, double, double>::value));
  EXPECT_TRUE((require_all_arithmetic_tester<int, int, int>::value));
}

/**
 * Require Not all arithmetic
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_all_not_arithmetic_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_all_not_arithmetic_tester<
    T1, T2, T3, stan::require_all_not_arithmetic<T1, T2, T3>> : std::true_type {
};

TEST(requires, not_all_arithmetic_test) {
  EXPECT_TRUE(
      (require_all_not_arithmetic_tester<double, std::string, double>::value));
  EXPECT_FALSE(
      (require_all_not_arithmetic_tester<double, double, double>::value));
  EXPECT_FALSE((require_all_not_arithmetic_tester<int, int, int>::value));
}

/**
 * Require any arithmetic
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_any_arithmetic_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_any_arithmetic_tester<T1, T2, T3,
                                     stan::require_any_arithmetic<T1, T2, T3>>
    : std::true_type {};

TEST(requires, any_arithmetic_test) {
  EXPECT_FALSE((require_any_arithmetic_tester<std::string, std::string,
                                              std::string>::value));
  EXPECT_TRUE(
      (require_any_arithmetic_tester<double, std::string, double>::value));
  EXPECT_TRUE((require_any_arithmetic_tester<int, int, int>::value));
}

/**
 * Require Not any arithmetic
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_any_not_arithmetic_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_any_not_arithmetic_tester<
    T1, T2, T3, stan::require_any_not_arithmetic<T1, T2, T3>> : std::true_type {
};

TEST(requires, not_any_arithmetic_test) {
  EXPECT_TRUE((require_any_not_arithmetic_tester<std::string, std::string,
                                                 std::string>::value));
  EXPECT_FALSE(
      (require_any_not_arithmetic_tester<double, double, double>::value));
  EXPECT_FALSE((require_any_not_arithmetic_tester<int, int, int>::value));
}

//  Require var_or_arithmetic

TEST(requires, var_or_arithmetic_test) {
  using stan::require_var_or_arithmetic;
  EXPECT_FALSE((unary_require_tester<require_var_or_arithmetic, std::string>::value));
  EXPECT_TRUE((unary_require_tester<require_var_or_arithmetic, double, int>::value));
}

// Require not var or arithmetic

TEST(requires, not_var_or_arithmetic_test) {
  using stan::require_not_var_or_arithmetic;
  EXPECT_TRUE((unary_require_tester<require_not_var_or_arithmetic, std::string>::value));
  EXPECT_FALSE((unary_require_tester<require_not_var_or_arithmetic, double, int>::value));
}

/**
 * Require all var_or_arithmetic
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_all_var_or_arithmetic_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_all_var_or_arithmetic_tester<
    T1, T2, T3, stan::require_all_var_or_arithmetic<T1, T2, T3>>
    : std::true_type {};

TEST(requires, all_var_or_arithmetic_test) {
  EXPECT_FALSE((require_all_var_or_arithmetic_tester<double, std::string,
                                                     double>::value));
  EXPECT_TRUE(
      (require_all_var_or_arithmetic_tester<double, double, double>::value));
  EXPECT_TRUE((require_all_var_or_arithmetic_tester<int, int, int>::value));
}

/**
 * Require Not all var_or_arithmetic
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_all_not_var_or_arithmetic_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_all_not_var_or_arithmetic_tester<
    T1, T2, T3, stan::require_all_not_var_or_arithmetic<T1, T2, T3>>
    : std::true_type {};

TEST(requires, not_all_var_or_arithmetic_test) {
  EXPECT_TRUE((require_all_not_var_or_arithmetic_tester<double, std::string,
                                                        double>::value));
  EXPECT_FALSE((
      require_all_not_var_or_arithmetic_tester<double, double, double>::value));
  EXPECT_FALSE(
      (require_all_not_var_or_arithmetic_tester<int, int, int>::value));
}

/**
 * Require any var_or_arithmetic
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_any_var_or_arithmetic_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_any_var_or_arithmetic_tester<
    T1, T2, T3, stan::require_any_var_or_arithmetic<T1, T2, T3>>
    : std::true_type {};

TEST(requires, any_var_or_arithmetic_test) {
  EXPECT_FALSE((require_any_var_or_arithmetic_tester<std::string, std::string,
                                                     std::string>::value));
  EXPECT_TRUE((require_any_var_or_arithmetic_tester<double, std::string,
                                                    double>::value));
  EXPECT_TRUE((require_any_var_or_arithmetic_tester<int, int, int>::value));
}

/**
 * Require Not any var_or_arithmetic
 */
template <typename T1, typename T2, typename T3, typename = void>
struct require_any_not_var_or_arithmetic_tester : std::false_type {};

template <typename T1, typename T2, typename T3>
struct require_any_not_var_or_arithmetic_tester<
    T1, T2, T3, stan::require_any_not_var_or_arithmetic<T1, T2, T3>>
    : std::true_type {};

TEST(requires, not_any_var_or_arithmetic_test) {
  EXPECT_TRUE(
      (require_any_not_var_or_arithmetic_tester<std::string, std::string,
                                                std::string>::value));
  EXPECT_FALSE((
      require_any_not_var_or_arithmetic_tester<double, double, double>::value));
  EXPECT_FALSE(
      (require_any_not_var_or_arithmetic_tester<int, int, int>::value));
}
