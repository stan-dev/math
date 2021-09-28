#ifdef STAN_OPENCL

#include <stan/math/opencl/prim.hpp>
#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <string>

/**
 * Require matrix_cl
 */
template <template <class> class TypeCheck, class Check, typename = void>
struct require_matrix_cl_tester : std::false_type {};

template <template <class> class TypeCheck, class Check>
struct require_matrix_cl_tester<TypeCheck, Check,
                                stan::require_matrix_cl_vt<TypeCheck, Check>>
    : std::true_type {};

TEST(requires, matrix_cl_test) {
  using stan::math::matrix_cl;
  EXPECT_FALSE(
      (require_matrix_cl_tester<std::is_floating_point, double>::value));
  EXPECT_TRUE(
      (require_matrix_cl_tester<std::is_arithmetic, matrix_cl<double>>::value));
  EXPECT_TRUE((require_matrix_cl_tester<std::is_floating_point,
                                        matrix_cl<double>>::value));
  EXPECT_FALSE((
      require_matrix_cl_tester<std::is_floating_point, matrix_cl<int>>::value));
}

/**
 * Require not matrix_cl
 */
template <template <class> class TypeCheck, class Check, typename = void>
struct require_not_matrix_cl_tester : std::false_type {};

template <template <class> class TypeCheck, class Check>
struct require_not_matrix_cl_tester<
    TypeCheck, Check, stan::require_not_matrix_cl_vt<TypeCheck, Check>>
    : std::true_type {};

TEST(requires, not_matrix_cl_test) {
  using stan::math::matrix_cl;
  EXPECT_TRUE(
      (require_not_matrix_cl_tester<std::is_floating_point, double>::value));
  EXPECT_FALSE((require_not_matrix_cl_tester<std::is_arithmetic,
                                             matrix_cl<double>>::value));
  EXPECT_FALSE((require_not_matrix_cl_tester<std::is_floating_point,
                                             matrix_cl<double>>::value));
  EXPECT_TRUE((require_not_matrix_cl_tester<std::is_floating_point,
                                            matrix_cl<int>>::value));
}

///////

/**
 * Require all matrix_cl
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_any_matrix_cl_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_any_matrix_cl_tester<
    TypeCheck, Check1, Check2,
    stan::require_any_matrix_cl_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, any_matrix_cl_test) {
  using stan::math::matrix_cl;
  EXPECT_TRUE(
      (require_any_matrix_cl_tester<std::is_floating_point, matrix_cl<double>,
                                    matrix_cl<double>>::value));
  EXPECT_TRUE((require_any_matrix_cl_tester<std::is_floating_point, double,
                                            matrix_cl<double>>::value));
  EXPECT_TRUE((require_any_matrix_cl_tester<std::is_floating_point,
                                            matrix_cl<double>, double>::value));
  EXPECT_TRUE(
      (require_any_matrix_cl_tester<std::is_floating_point, matrix_cl<int>,
                                    matrix_cl<double>>::value));
  EXPECT_FALSE((require_any_matrix_cl_tester<std::is_floating_point, int,
                                             std::string>::value));
  EXPECT_FALSE((
      require_any_matrix_cl_tester<std::is_arithmetic, double, double>::value));
}

////////

/**
 * Require all matrix_cl
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_any_not_matrix_cl_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_any_not_matrix_cl_tester<
    TypeCheck, Check1, Check2,
    stan::require_any_not_matrix_cl_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, any_not_matrix_cl_test) {
  using stan::math::matrix_cl;
  EXPECT_FALSE((require_any_not_matrix_cl_tester<std::is_floating_point,
                                                 matrix_cl<double>,
                                                 matrix_cl<double>>::value));
  EXPECT_TRUE((require_any_not_matrix_cl_tester<std::is_floating_point, double,
                                                matrix_cl<double>>::value));
  EXPECT_TRUE(
      (require_any_not_matrix_cl_tester<std::is_floating_point,
                                        matrix_cl<double>, double>::value));
  EXPECT_TRUE(
      (require_any_not_matrix_cl_tester<std::is_floating_point, matrix_cl<int>,
                                        matrix_cl<double>>::value));
  EXPECT_TRUE((require_any_not_matrix_cl_tester<std::is_floating_point, int,
                                                std::string>::value));
  EXPECT_TRUE((require_any_not_matrix_cl_tester<std::is_arithmetic, double,
                                                double>::value));
}

/////

/**
 * Require all matrix_cl
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_all_matrix_cl_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_all_matrix_cl_tester<
    TypeCheck, Check1, Check2,
    stan::require_all_matrix_cl_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, all_matrix_cl_test) {
  using stan::math::matrix_cl;
  EXPECT_TRUE(
      (require_all_matrix_cl_tester<std::is_floating_point, matrix_cl<double>,
                                    matrix_cl<double>>::value));
  EXPECT_FALSE((require_all_matrix_cl_tester<std::is_floating_point, double,
                                             matrix_cl<double>>::value));
  EXPECT_FALSE(
      (require_all_matrix_cl_tester<std::is_floating_point, matrix_cl<double>,
                                    double>::value));
  EXPECT_FALSE(
      (require_all_matrix_cl_tester<std::is_floating_point, matrix_cl<int>,
                                    matrix_cl<double>>::value));
  EXPECT_FALSE((require_all_matrix_cl_tester<std::is_floating_point, int,
                                             std::string>::value));
  EXPECT_FALSE((
      require_all_matrix_cl_tester<std::is_arithmetic, double, double>::value));
}

/**
 * Require all not matrix_cl
 */
template <template <class> class TypeCheck, class Check1, class Check2,
          typename = void>
struct require_all_not_matrix_cl_tester : std::false_type {};

template <template <class> class TypeCheck, class Check1, class Check2>
struct require_all_not_matrix_cl_tester<
    TypeCheck, Check1, Check2,
    stan::require_all_not_matrix_cl_vt<TypeCheck, Check1, Check2>>
    : std::true_type {};

TEST(requires, all_not_matrix_cl_test) {
  using stan::math::matrix_cl;
  EXPECT_FALSE((require_all_not_matrix_cl_tester<std::is_floating_point,
                                                 matrix_cl<double>,
                                                 matrix_cl<double>>::value));
  EXPECT_FALSE((require_all_not_matrix_cl_tester<std::is_floating_point, double,
                                                 matrix_cl<double>>::value));
  EXPECT_FALSE(
      (require_all_not_matrix_cl_tester<std::is_floating_point,
                                        matrix_cl<double>, double>::value));
  EXPECT_FALSE(
      (require_all_not_matrix_cl_tester<std::is_floating_point, matrix_cl<int>,
                                        matrix_cl<double>>::value));
  EXPECT_TRUE((require_all_not_matrix_cl_tester<std::is_floating_point, int,
                                                std::string>::value));
  EXPECT_TRUE((require_all_not_matrix_cl_tester<std::is_arithmetic, double,
                                                double>::value));
}
#endif  // STAN_OPENCL
