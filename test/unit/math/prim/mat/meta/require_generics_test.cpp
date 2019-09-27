#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <string>

/**
 * Require eigen
 */
template <template <class> class TypeCheck, class Check, typename = void>
struct require_eigen_tester : std::false_type {};

template <template <class> class TypeCheck, class Check>
struct require_eigen_tester<TypeCheck, Check,
                            stan::require_eigen_vt<TypeCheck, Check>>
    : std::true_type {};

TEST(requires, eigen_test) {
  EXPECT_FALSE((require_eigen_tester<std::is_floating_point, double>::value));
  EXPECT_TRUE((require_eigen_tester<std::is_arithmetic,
                                    Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_TRUE((require_eigen_tester<std::is_floating_point,
                                    Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_FALSE(
      (require_eigen_tester<std::is_floating_point,
                            Eigen::Matrix<std::string, -1, -1>>::value));
}

/**
 * Require not eigen
 */
template <template <class> class TypeCheck, class Check, typename = void>
struct require_not_eigen_tester : std::false_type {};

template <template <class> class TypeCheck, class Check>
struct require_not_eigen_tester<TypeCheck, Check,
                                stan::require_not_eigen_vt<TypeCheck, Check>>
    : std::true_type {};

TEST(requires, not_eigen_test) {
  EXPECT_TRUE(
      (require_not_eigen_tester<std::is_floating_point, double>::value));
  EXPECT_FALSE(
      (require_not_eigen_tester<std::is_arithmetic,
                                Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_FALSE(
      (require_not_eigen_tester<std::is_floating_point,
                                Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_TRUE(
      (require_not_eigen_tester<std::is_floating_point,
                                Eigen::Matrix<std::string, -1, -1>>::value));
}

///////

/**
 * Require all eigen
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_any_eigen_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_any_eigen_tester<
    TypeCheck, Check1, Check2,
    stan::require_any_eigen_vt<TypeCheck, Check1, Check2>> : std::true_type {};

TEST(requires, any_eigen_test) {
  EXPECT_TRUE((require_any_eigen_tester<std::is_floating_point,
                                        Eigen::Matrix<double, -1, -1>,
                                        Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_TRUE((require_any_eigen_tester<std::is_floating_point, double,
                                        Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_TRUE(
      (require_any_eigen_tester<std::is_floating_point,
                                Eigen::Matrix<double, -1, -1>, double>::value));
  EXPECT_TRUE((require_any_eigen_tester<std::is_floating_point,
                                        Eigen::Matrix<std::string, -1, -1>,
                                        Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_FALSE((require_any_eigen_tester<std::is_floating_point, int,
                                         std::string>::value));
  EXPECT_FALSE(
      (require_any_eigen_tester<std::is_arithmetic, double, double>::value));
}

////////

/**
 * Require all eigen
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_any_not_eigen_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_any_not_eigen_tester<
    TypeCheck, Check1, Check2,
    stan::require_any_not_eigen_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, not_any_eigen_test) {
  EXPECT_FALSE(
      (require_any_not_eigen_tester<std::is_floating_point,
                                    Eigen::Matrix<double, -1, -1>,
                                    Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_FALSE(
      (require_any_not_eigen_tester<std::is_floating_point, double,
                                    Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_FALSE((require_any_not_eigen_tester<std::is_floating_point,
                                             Eigen::Matrix<double, -1, -1>,
                                             double>::value));
  EXPECT_FALSE(
      (require_any_not_eigen_tester<std::is_floating_point,
                                    Eigen::Matrix<std::string, -1, -1>,
                                    Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_TRUE((require_any_not_eigen_tester<std::is_floating_point, int,
                                            std::string>::value));
  EXPECT_TRUE((
      require_any_not_eigen_tester<std::is_arithmetic, double, double>::value));
}

/////

/**
 * Require all eigen
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_all_eigen_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_all_eigen_tester<
    TypeCheck, Check1, Check2,
    stan::require_all_eigen_vt<TypeCheck, Check1, Check2>> : std::true_type {};

TEST(requires, all_eigen_test) {
  EXPECT_TRUE((require_all_eigen_tester<std::is_floating_point,
                                        Eigen::Matrix<double, -1, -1>,
                                        Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_FALSE(
      (require_all_eigen_tester<std::is_floating_point, double,
                                Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_FALSE(
      (require_all_eigen_tester<std::is_floating_point,
                                Eigen::Matrix<double, -1, -1>, double>::value));
  EXPECT_FALSE(
      (require_all_eigen_tester<std::is_floating_point,
                                Eigen::Matrix<std::string, -1, -1>,
                                Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_FALSE((require_all_eigen_tester<std::is_floating_point, int,
                                         std::string>::value));
  EXPECT_FALSE(
      (require_all_eigen_tester<std::is_arithmetic, double, double>::value));
}

/**
 * Require all not eigen
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_all_not_eigen_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_all_not_eigen_tester<
    TypeCheck, Check1, Check2,
    stan::require_all_not_eigen_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, not_all_eigen_test) {
  EXPECT_FALSE(
      (require_all_not_eigen_tester<std::is_floating_point,
                                    Eigen::Matrix<double, -1, -1>,
                                    Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_TRUE(
      (require_all_not_eigen_tester<std::is_floating_point, double,
                                    Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_TRUE((require_all_not_eigen_tester<std::is_floating_point,
                                            Eigen::Matrix<double, -1, -1>,
                                            double>::value));
  EXPECT_TRUE(
      (require_all_not_eigen_tester<std::is_floating_point,
                                    Eigen::Matrix<std::string, -1, -1>,
                                    Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_TRUE((require_all_not_eigen_tester<std::is_floating_point, int,
                                            std::string>::value));
  EXPECT_TRUE((
      require_all_not_eigen_tester<std::is_arithmetic, double, double>::value));
}

/**
 * Require eigen vector
 */
template <template <class> class TypeCheck, class Check, typename = void>
struct require_eigen_vector_tester : std::false_type {};

template <template <class> class TypeCheck, class Check>
struct require_eigen_vector_tester<
    TypeCheck, Check, stan::require_eigen_vector_vt<TypeCheck, Check>>
    : std::true_type {};

TEST(requires, eigen_vector_test) {
  EXPECT_FALSE(
      (require_eigen_vector_tester<std::is_floating_point, double>::value));
  EXPECT_TRUE(
      (require_eigen_vector_tester<std::is_arithmetic,
                                   Eigen::Matrix<double, 1, -1>>::value));
  EXPECT_TRUE(
      (require_eigen_vector_tester<std::is_floating_point,
                                   Eigen::Matrix<double, 1, -1>>::value));
  EXPECT_FALSE(
      (require_eigen_vector_tester<std::is_floating_point,
                                   Eigen::Matrix<std::string, 1, -1>>::value));
}

/**
 * Require not eigen
 */
template <template <class> class TypeCheck, class Check, typename = void>
struct require_not_eigen_vector_tester : std::false_type {};

template <template <class> class TypeCheck, class Check>
struct require_not_eigen_vector_tester<
    TypeCheck, Check, stan::require_not_eigen_vector_vt<TypeCheck, Check>>
    : std::true_type {};

TEST(requires, not_eigen_vector_test) {
  EXPECT_TRUE(
      (require_not_eigen_vector_tester<std::is_floating_point, double>::value));
  EXPECT_FALSE(
      (require_not_eigen_vector_tester<std::is_arithmetic,
                                       Eigen::Matrix<double, 1, -1>>::value));
  EXPECT_FALSE(
      (require_not_eigen_vector_tester<std::is_floating_point,
                                       Eigen::Matrix<double, 1, -1>>::value));
  EXPECT_TRUE(
      (require_not_eigen_vector_tester<
          std::is_floating_point, Eigen::Matrix<std::string, 1, -1>>::value));
}

///////

/**
 * Require all eigen
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_any_eigen_vector_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_any_eigen_vector_tester<
    TypeCheck, Check1, Check2,
    stan::require_any_eigen_vector_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, any_eigen_vector_test) {
  EXPECT_TRUE(
      (require_any_eigen_vector_tester<std::is_floating_point,
                                       Eigen::Matrix<double, 1, -1>,
                                       Eigen::Matrix<double, 1, -1>>::value));
  EXPECT_TRUE(
      (require_any_eigen_vector_tester<std::is_floating_point, double,
                                       Eigen::Matrix<double, 1, -1>>::value));
  EXPECT_TRUE((require_any_eigen_vector_tester<std::is_floating_point,
                                               Eigen::Matrix<double, 1, -1>,
                                               double>::value));
  EXPECT_TRUE(
      (require_any_eigen_vector_tester<std::is_floating_point,
                                       Eigen::Matrix<std::string, 1, -1>,
                                       Eigen::Matrix<double, 1, -1>>::value));
  EXPECT_FALSE((require_any_eigen_vector_tester<std::is_floating_point, int,
                                                std::string>::value));
  EXPECT_FALSE((require_any_eigen_vector_tester<std::is_arithmetic, double,
                                                double>::value));
}

////////

/**
 * Require all eigen
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_any_not_eigen_vector_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_any_not_eigen_vector_tester<
    TypeCheck, Check1, Check2,
    stan::require_any_not_eigen_vector_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, not_any_eigen_vector_test) {
  EXPECT_FALSE((require_any_not_eigen_vector_tester<
                std::is_floating_point, Eigen::Matrix<double, 1, -1>,
                Eigen::Matrix<double, 1, -1>>::value));
  EXPECT_FALSE((require_any_not_eigen_vector_tester<
                std::is_floating_point, double,
                Eigen::Matrix<double, 1, -1>>::value));
  EXPECT_FALSE(
      (require_any_not_eigen_vector_tester<std::is_floating_point,
                                           Eigen::Matrix<double, 1, -1>,
                                           double>::value));
  EXPECT_FALSE((require_any_not_eigen_vector_tester<
                std::is_floating_point, Eigen::Matrix<std::string, 1, -1>,
                Eigen::Matrix<double, 1, -1>>::value));
  EXPECT_TRUE((require_any_not_eigen_vector_tester<std::is_floating_point, int,
                                                   std::string>::value));
  EXPECT_TRUE((require_any_not_eigen_vector_tester<std::is_arithmetic, double,
                                                   double>::value));
}

/////

/**
 * Require all eigen
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_all_eigen_vector_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_all_eigen_vector_tester<
    TypeCheck, Check1, Check2,
    stan::require_all_eigen_vector_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, all_eigen_vector_test) {
  EXPECT_TRUE(
      (require_all_eigen_vector_tester<std::is_floating_point,
                                       Eigen::Matrix<double, 1, -1>,
                                       Eigen::Matrix<double, 1, -1>>::value));
  EXPECT_FALSE(
      (require_all_eigen_vector_tester<std::is_floating_point, double,
                                       Eigen::Matrix<double, 1, -1>>::value));
  EXPECT_FALSE((require_all_eigen_vector_tester<std::is_floating_point,
                                                Eigen::Matrix<double, 1, -1>,
                                                double>::value));
  EXPECT_FALSE(
      (require_all_eigen_vector_tester<std::is_floating_point,
                                       Eigen::Matrix<std::string, 1, -1>,
                                       Eigen::Matrix<double, 1, -1>>::value));
  EXPECT_FALSE((require_all_eigen_vector_tester<std::is_floating_point, int,
                                                std::string>::value));
  EXPECT_FALSE((require_all_eigen_vector_tester<std::is_arithmetic, double,
                                                double>::value));
}

/**
 * Require all not eigen
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_all_not_eigen_vector_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_all_not_eigen_vector_tester<
    TypeCheck, Check1, Check2,
    stan::require_all_not_eigen_vector_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, not_all_eigen_vector_test) {
  EXPECT_FALSE((require_all_not_eigen_vector_tester<
                std::is_floating_point, Eigen::Matrix<double, 1, -1>,
                Eigen::Matrix<double, 1, -1>>::value));
  EXPECT_TRUE((require_all_not_eigen_vector_tester<
               std::is_floating_point, double,
               Eigen::Matrix<double, 1, -1>>::value));
  EXPECT_TRUE((require_all_not_eigen_vector_tester<std::is_floating_point,
                                                   Eigen::Matrix<double, 1, -1>,
                                                   double>::value));
  EXPECT_TRUE((require_all_not_eigen_vector_tester<
               std::is_floating_point, Eigen::Matrix<std::string, 1, -1>,
               Eigen::Matrix<double, 1, -1>>::value));
  EXPECT_TRUE((require_all_not_eigen_vector_tester<std::is_floating_point, int,
                                                   std::string>::value));
  EXPECT_TRUE((require_all_not_eigen_vector_tester<std::is_arithmetic, double,
                                                   double>::value));
}
