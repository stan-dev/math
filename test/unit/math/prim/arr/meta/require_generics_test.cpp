#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <vector>
#include <string>

TEST(requires, container_type_test) {
  EXPECT_FALSE((stan::container_value_type_check_base<
                stan::is_vector, std::is_floating_point, std::string>::value));
  EXPECT_TRUE(
      (stan::container_value_type_check_base<stan::is_vector,
                                             std::is_floating_point,
                                             std::vector<double>>::value));
}

////////////////////////////////
/**
 * Require container
 */
template <template <class> class ContainerCheck,
          template <class> class TypeCheck, class Check, typename = void>
struct require_container_tester : std::false_type {};

template <template <class> class ContainerCheck,
          template <class> class TypeCheck, class Check>
struct require_container_tester<
    ContainerCheck, TypeCheck, Check,
    stan::require_container_vt<ContainerCheck, TypeCheck, Check>>
    : std::true_type {};

TEST(requires, generic_container_type_test) {
  EXPECT_FALSE(
      (require_container_tester<stan::is_vector, std::is_floating_point,
                                double>::value));
  EXPECT_TRUE((require_container_tester<stan::is_vector, std::is_floating_point,
                                        std::vector<double>>::value));
}

/**
 * Require std_vector
 */
template <template <class> class TypeCheck, class Check, typename = void>
struct require_std_vector_tester : std::false_type {};

template <template <class> class TypeCheck, class Check>
struct require_std_vector_tester<TypeCheck, Check,
                                 stan::require_std_vector_vt<TypeCheck, Check>>
    : std::true_type {};

TEST(requires, std_vector_test) {
  EXPECT_FALSE(
      (require_std_vector_tester<std::is_floating_point, double>::value));
  EXPECT_TRUE((require_std_vector_tester<std::is_arithmetic,
                                         std::vector<double>>::value));
  EXPECT_TRUE((require_std_vector_tester<std::is_floating_point,
                                         std::vector<double>>::value));
  EXPECT_FALSE((require_std_vector_tester<std::is_floating_point,
                                          std::vector<std::string>>::value));
}

/**
 * Require not std_vector
 */
template <template <class> class TypeCheck, class Check, typename = void>
struct require_not_std_vector_tester : std::false_type {};

template <template <class> class TypeCheck, class Check>
struct require_not_std_vector_tester<
    TypeCheck, Check, stan::require_not_std_vector_vt<TypeCheck, Check>>
    : std::true_type {};

TEST(requires, not_std_vector_test) {
  EXPECT_TRUE(
      (require_not_std_vector_tester<std::is_floating_point, double>::value));
  EXPECT_FALSE((require_not_std_vector_tester<std::is_arithmetic,
                                              std::vector<double>>::value));
  EXPECT_FALSE((require_not_std_vector_tester<std::is_floating_point,
                                              std::vector<double>>::value));
  EXPECT_TRUE((require_not_std_vector_tester<std::is_floating_point,
                                             std::vector<std::string>>::value));
}

///////

/**
 * Require all std_vector
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_any_std_vector_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_any_std_vector_tester<
    TypeCheck, Check1, Check2,
    stan::require_any_std_vector_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, any_std_vector_test) {
  EXPECT_TRUE((
      require_any_std_vector_tester<std::is_floating_point, std::vector<double>,
                                    std::vector<double>>::value));
  EXPECT_TRUE((require_any_std_vector_tester<std::is_floating_point, double,
                                             std::vector<double>>::value));
  EXPECT_TRUE(
      (require_any_std_vector_tester<std::is_floating_point,
                                     std::vector<double>, double>::value));
  EXPECT_TRUE((require_any_std_vector_tester<std::is_floating_point,
                                             std::vector<std::string>,
                                             std::vector<double>>::value));
  EXPECT_FALSE((require_any_std_vector_tester<std::is_floating_point, int,
                                              std::string>::value));
  EXPECT_FALSE((require_any_std_vector_tester<std::is_arithmetic, double,
                                              double>::value));
}

////////

/**
 * Require all std_vector
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_any_not_std_vector_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_any_not_std_vector_tester<
    TypeCheck, Check1, Check2,
    stan::require_any_not_std_vector_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, not_any_std_vector_test) {
  EXPECT_FALSE((require_any_not_std_vector_tester<std::is_floating_point,
                                                  std::vector<double>,
                                                  std::vector<double>>::value));
  EXPECT_FALSE(
      (require_any_not_std_vector_tester<std::is_floating_point, double,
                                         std::vector<double>>::value));
  EXPECT_FALSE(
      (require_any_not_std_vector_tester<std::is_floating_point,
                                         std::vector<double>, double>::value));
  EXPECT_FALSE((require_any_not_std_vector_tester<std::is_floating_point,
                                                  std::vector<std::string>,
                                                  std::vector<double>>::value));
  EXPECT_TRUE((require_any_not_std_vector_tester<std::is_floating_point, int,
                                                 std::string>::value));
  EXPECT_TRUE((require_any_not_std_vector_tester<std::is_arithmetic, double,
                                                 double>::value));
}

/////

/**
 * Require all std_vector
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_all_std_vector_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_all_std_vector_tester<
    TypeCheck, Check1, Check2,
    stan::require_all_std_vector_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, all_std_vector_test) {
  EXPECT_TRUE((
      require_all_std_vector_tester<std::is_floating_point, std::vector<double>,
                                    std::vector<double>>::value));
  EXPECT_FALSE((require_all_std_vector_tester<std::is_floating_point, double,
                                              std::vector<double>>::value));
  EXPECT_FALSE(
      (require_all_std_vector_tester<std::is_floating_point,
                                     std::vector<double>, double>::value));
  EXPECT_FALSE((require_all_std_vector_tester<std::is_floating_point,
                                              std::vector<std::string>,
                                              std::vector<double>>::value));
  EXPECT_FALSE((require_all_std_vector_tester<std::is_floating_point, int,
                                              std::string>::value));
  EXPECT_FALSE((require_all_std_vector_tester<std::is_arithmetic, double,
                                              double>::value));
}

/**
 * Require all not std_vector
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_all_not_std_vector_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_all_not_std_vector_tester<
    TypeCheck, Check1, Check2,
    stan::require_all_not_std_vector_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, not_all_std_vector_test) {
  EXPECT_FALSE((require_all_not_std_vector_tester<std::is_floating_point,
                                                  std::vector<double>,
                                                  std::vector<double>>::value));
  EXPECT_TRUE((require_all_not_std_vector_tester<std::is_floating_point, double,
                                                 std::vector<double>>::value));
  EXPECT_TRUE(
      (require_all_not_std_vector_tester<std::is_floating_point,
                                         std::vector<double>, double>::value));
  EXPECT_TRUE((require_all_not_std_vector_tester<std::is_floating_point,
                                                 std::vector<std::string>,
                                                 std::vector<double>>::value));
  EXPECT_TRUE((require_all_not_std_vector_tester<std::is_floating_point, int,
                                                 std::string>::value));
  EXPECT_TRUE((require_all_not_std_vector_tester<std::is_arithmetic, double,
                                                 double>::value));
}

/**
 * Require vector
 */
template <template <class> class TypeCheck, class Check, typename = void>
struct require_vector_tester : std::false_type {};

template <template <class> class TypeCheck, class Check>
struct require_vector_tester<TypeCheck, Check,
                             stan::require_vector_vt<TypeCheck, Check>>
    : std::true_type {};

TEST(requires, vector_test) {
  EXPECT_FALSE((require_vector_tester<std::is_floating_point, double>::value));
  EXPECT_TRUE(
      (require_vector_tester<std::is_arithmetic, std::vector<double>>::value));
  EXPECT_TRUE((require_vector_tester<std::is_floating_point,
                                     std::vector<double>>::value));
  EXPECT_FALSE((require_vector_tester<std::is_floating_point,
                                      std::vector<std::string>>::value));
}

/**
 * Require not vector
 */
template <template <class> class TypeCheck, class Check, typename = void>
struct require_not_vector_tester : std::false_type {};

template <template <class> class TypeCheck, class Check>
struct require_not_vector_tester<TypeCheck, Check,
                                 stan::require_not_vector_vt<TypeCheck, Check>>
    : std::true_type {};

TEST(requires, not_vector_test) {
  EXPECT_TRUE(
      (require_not_vector_tester<std::is_floating_point, double>::value));
  EXPECT_FALSE((require_not_vector_tester<std::is_arithmetic,
                                          std::vector<double>>::value));
  EXPECT_FALSE((require_not_vector_tester<std::is_floating_point,
                                          std::vector<double>>::value));
  EXPECT_TRUE((require_not_vector_tester<std::is_floating_point,
                                         std::vector<std::string>>::value));
}

///////

/**
 * Require all vector
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_any_vector_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_any_vector_tester<
    TypeCheck, Check1, Check2,
    stan::require_any_vector_vt<TypeCheck, Check1, Check2>> : std::true_type {};

TEST(requires, any_vector_test) {
  EXPECT_TRUE(
      (require_any_vector_tester<std::is_floating_point, std::vector<double>,
                                 std::vector<double>>::value));
  EXPECT_TRUE((require_any_vector_tester<std::is_floating_point, double,
                                         std::vector<double>>::value));
  EXPECT_TRUE((require_any_vector_tester<std::is_floating_point,
                                         std::vector<double>, double>::value));
  EXPECT_TRUE((require_any_vector_tester<std::is_floating_point,
                                         std::vector<std::string>,
                                         std::vector<double>>::value));
  EXPECT_FALSE((require_any_vector_tester<std::is_floating_point, int,
                                          std::string>::value));
  EXPECT_FALSE(
      (require_any_vector_tester<std::is_arithmetic, double, double>::value));
}

////////

/**
 * Require all vector
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_any_not_vector_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_any_not_vector_tester<
    TypeCheck, Check1, Check2,
    stan::require_any_not_vector_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, not_any_vector_test) {
  EXPECT_FALSE((
      require_any_not_vector_tester<std::is_floating_point, std::vector<double>,
                                    std::vector<double>>::value));
  EXPECT_FALSE((require_any_not_vector_tester<std::is_floating_point, double,
                                              std::vector<double>>::value));
  EXPECT_FALSE(
      (require_any_not_vector_tester<std::is_floating_point,
                                     std::vector<double>, double>::value));
  EXPECT_FALSE((require_any_not_vector_tester<std::is_floating_point,
                                              std::vector<std::string>,
                                              std::vector<double>>::value));
  EXPECT_TRUE((require_any_not_vector_tester<std::is_floating_point, int,
                                             std::string>::value));
  EXPECT_TRUE((require_any_not_vector_tester<std::is_arithmetic, double,
                                             double>::value));
}

/////

/**
 * Require all vector
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_all_vector_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_all_vector_tester<
    TypeCheck, Check1, Check2,
    stan::require_all_vector_vt<TypeCheck, Check1, Check2>> : std::true_type {};

TEST(requires, all_vector_test) {
  EXPECT_TRUE(
      (require_all_vector_tester<std::is_floating_point, std::vector<double>,
                                 std::vector<double>>::value));
  EXPECT_FALSE((require_all_vector_tester<std::is_floating_point, double,
                                          std::vector<double>>::value));
  EXPECT_FALSE((require_all_vector_tester<std::is_floating_point,
                                          std::vector<double>, double>::value));
  EXPECT_FALSE((require_all_vector_tester<std::is_floating_point,
                                          std::vector<std::string>,
                                          std::vector<double>>::value));
  EXPECT_FALSE((require_all_vector_tester<std::is_floating_point, int,
                                          std::string>::value));
  EXPECT_FALSE(
      (require_all_vector_tester<std::is_arithmetic, double, double>::value));
}

/**
 * Require all not vector
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_all_not_vector_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_all_not_vector_tester<
    TypeCheck, Check1, Check2,
    stan::require_all_not_vector_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, not_all_vector_test) {
  EXPECT_FALSE((
      require_all_not_vector_tester<std::is_floating_point, std::vector<double>,
                                    std::vector<double>>::value));
  EXPECT_TRUE((require_all_not_vector_tester<std::is_floating_point, double,
                                             std::vector<double>>::value));
  EXPECT_TRUE(
      (require_all_not_vector_tester<std::is_floating_point,
                                     std::vector<double>, double>::value));
  EXPECT_TRUE((require_all_not_vector_tester<std::is_floating_point,
                                             std::vector<std::string>,
                                             std::vector<double>>::value));
  EXPECT_TRUE((require_all_not_vector_tester<std::is_floating_point, int,
                                             std::string>::value));
  EXPECT_TRUE((require_all_not_vector_tester<std::is_arithmetic, double,
                                             double>::value));
}

/**
 * Require vector
 */
template <template <class> class TypeCheck, class Check, typename = void>
struct require_vector_like_tester : std::false_type {};

template <template <class> class TypeCheck, class Check>
struct require_vector_like_tester<
    TypeCheck, Check, stan::require_vector_like_vt<TypeCheck, Check>>
    : std::true_type {};

TEST(requires, vector_like_test) {
  EXPECT_FALSE(
      (require_vector_like_tester<std::is_floating_point, double>::value));
  EXPECT_TRUE((require_vector_like_tester<std::is_arithmetic,
                                          std::vector<double>>::value));
  EXPECT_TRUE((require_vector_like_tester<std::is_floating_point,
                                          std::vector<double>>::value));
  EXPECT_FALSE((require_vector_like_tester<std::is_floating_point,
                                           std::vector<std::string>>::value));
}

/**
 * Require not vector
 */
template <template <class> class TypeCheck, class Check, typename = void>
struct require_not_vector_like_tester : std::false_type {};

template <template <class> class TypeCheck, class Check>
struct require_not_vector_like_tester<
    TypeCheck, Check, stan::require_not_vector_like_vt<TypeCheck, Check>>
    : std::true_type {};

TEST(requires, not_vector_like_test) {
  EXPECT_TRUE(
      (require_not_vector_like_tester<std::is_floating_point, double>::value));
  EXPECT_FALSE((require_not_vector_like_tester<std::is_arithmetic,
                                               std::vector<double>>::value));
  EXPECT_FALSE((require_not_vector_like_tester<std::is_floating_point,
                                               std::vector<double>>::value));
  EXPECT_TRUE(
      (require_not_vector_like_tester<std::is_floating_point,
                                      std::vector<std::string>>::value));
}

///////

/**
 * Require all vector
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_any_vector_like_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_any_vector_like_tester<
    TypeCheck, Check1, Check2,
    stan::require_any_vector_like_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, any_vector_like_test) {
  EXPECT_TRUE((require_any_vector_like_tester<std::is_floating_point,
                                              std::vector<double>,
                                              std::vector<double>>::value));
  EXPECT_TRUE((require_any_vector_like_tester<std::is_floating_point, double,
                                              std::vector<double>>::value));
  EXPECT_TRUE(
      (require_any_vector_like_tester<std::is_floating_point,
                                      std::vector<double>, double>::value));
  EXPECT_TRUE((require_any_vector_like_tester<std::is_floating_point,
                                              std::vector<std::string>,
                                              std::vector<double>>::value));
  EXPECT_FALSE((require_any_vector_like_tester<std::is_floating_point, int,
                                               std::string>::value));
  EXPECT_FALSE((require_any_vector_like_tester<std::is_arithmetic, double,
                                               double>::value));
}

////////

/**
 * Require all vector
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_any_not_vector_like_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_any_not_vector_like_tester<
    TypeCheck, Check1, Check2,
    stan::require_any_not_vector_like_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, not_any_vector_like_test) {
  EXPECT_FALSE(
      (require_any_not_vector_like_tester<std::is_floating_point,
                                          std::vector<double>,
                                          std::vector<double>>::value));
  EXPECT_FALSE(
      (require_any_not_vector_like_tester<std::is_floating_point, double,
                                          std::vector<double>>::value));
  EXPECT_FALSE(
      (require_any_not_vector_like_tester<std::is_floating_point,
                                          std::vector<double>, double>::value));
  EXPECT_FALSE(
      (require_any_not_vector_like_tester<std::is_floating_point,
                                          std::vector<std::string>,
                                          std::vector<double>>::value));
  EXPECT_TRUE((require_any_not_vector_like_tester<std::is_floating_point, int,
                                                  std::string>::value));
  EXPECT_TRUE((require_any_not_vector_like_tester<std::is_arithmetic, double,
                                                  double>::value));
}

/////

/**
 * Require all vector
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_all_vector_like_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_all_vector_like_tester<
    TypeCheck, Check1, Check2,
    stan::require_all_vector_like_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, all_vector_like_test) {
  EXPECT_TRUE((require_all_vector_like_tester<std::is_floating_point,
                                              std::vector<double>,
                                              std::vector<double>>::value));
  EXPECT_FALSE((require_all_vector_like_tester<std::is_floating_point, double,
                                               std::vector<double>>::value));
  EXPECT_FALSE(
      (require_all_vector_like_tester<std::is_floating_point,
                                      std::vector<double>, double>::value));
  EXPECT_FALSE((require_all_vector_like_tester<std::is_floating_point,
                                               std::vector<std::string>,
                                               std::vector<double>>::value));
  EXPECT_FALSE((require_all_vector_like_tester<std::is_floating_point, int,
                                               std::string>::value));
  EXPECT_FALSE((require_all_vector_like_tester<std::is_arithmetic, double,
                                               double>::value));
}

/**
 * Require all not vector
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_all_not_vector_like_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_all_not_vector_like_tester<
    TypeCheck, Check1, Check2,
    stan::require_all_not_vector_like_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, not_all_vector_like_test) {
  EXPECT_FALSE(
      (require_all_not_vector_like_tester<std::is_floating_point,
                                          std::vector<double>,
                                          std::vector<double>>::value));
  EXPECT_TRUE(
      (require_all_not_vector_like_tester<std::is_floating_point, double,
                                          std::vector<double>>::value));
  EXPECT_TRUE(
      (require_all_not_vector_like_tester<std::is_floating_point,
                                          std::vector<double>, double>::value));
  EXPECT_TRUE((require_all_not_vector_like_tester<std::is_floating_point,
                                                  std::vector<std::string>,
                                                  std::vector<double>>::value));
  EXPECT_TRUE((require_all_not_vector_like_tester<std::is_floating_point, int,
                                                  std::string>::value));
  EXPECT_TRUE((require_all_not_vector_like_tester<std::is_arithmetic, double,
                                                  double>::value));
}
