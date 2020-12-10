#ifdef STAN_OPENCL

#include <stan/math/prim.hpp>
#include <stan/math/opencl/kernel_generator/elt_function_cl.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <test/unit/math/opencl/kernel_generator/reference_kernel.hpp>
#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <gtest/gtest.h>
#include <string>

using Eigen::Dynamic;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using stan::math::matrix_cl;

namespace stan {
namespace math {

MatrixXd rsqrt(const MatrixXd& a) { return stan::math::inv_sqrt(a); }

}  // namespace math
}  // namespace stan

#define TEST_UNARY_FUNCTION(fun)                       \
  TEST(KernelGenerator, fun##_test) {                  \
    using stan::math::fun;                             \
    MatrixXd m1(3, 3);                                 \
    m1 << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9; \
                                                       \
    matrix_cl<double> m1_cl(m1);                       \
    auto tmp = fun(m1_cl);                             \
    matrix_cl<double> res_cl = tmp;                    \
                                                       \
    MatrixXd res = stan::math::from_matrix_cl(res_cl); \
    MatrixXd correct = fun(m1);                        \
    EXPECT_NEAR_REL(correct, res);                     \
  }

TEST_UNARY_FUNCTION(rsqrt)
TEST_UNARY_FUNCTION(sqrt)
TEST_UNARY_FUNCTION(cbrt)

TEST_UNARY_FUNCTION(exp)
TEST_UNARY_FUNCTION(exp2)
TEST_UNARY_FUNCTION(expm1)
TEST_UNARY_FUNCTION(log)
TEST_UNARY_FUNCTION(log2)
TEST_UNARY_FUNCTION(log10)
TEST_UNARY_FUNCTION(log1p)

TEST_UNARY_FUNCTION(sin)
TEST_UNARY_FUNCTION(sinh)
TEST_UNARY_FUNCTION(cos)
TEST_UNARY_FUNCTION(cosh)
TEST_UNARY_FUNCTION(tan)
TEST_UNARY_FUNCTION(tanh)

TEST_UNARY_FUNCTION(asin)
TEST_UNARY_FUNCTION(asinh)
TEST_UNARY_FUNCTION(acos)

TEST(KernelGenerator, acosh_test) {
  std::string kernel_filename = "unary_function_acosh.cl";
  MatrixXd m1(3, 3);
  m1 << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9;

  matrix_cl<double> m1_cl(m1);
  auto tmp = stan::math::acosh(m1_cl);
  matrix_cl<double> res_cl;

  std::string kernel_src = tmp.get_kernel_source_for_evaluating_into(res_cl);
  stan::test::store_reference_kernel_if_needed(kernel_filename, kernel_src);
  std::string expected_kernel_src
      = stan::test::load_reference_kernel(kernel_filename);
  EXPECT_EQ(expected_kernel_src, kernel_src);

  res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = stan::math::acosh(m1);
  EXPECT_NEAR_REL(correct, res);
}

TEST_UNARY_FUNCTION(atan)
TEST_UNARY_FUNCTION(atanh)

TEST_UNARY_FUNCTION(tgamma)
TEST_UNARY_FUNCTION(lgamma)
TEST_UNARY_FUNCTION(erf)
TEST_UNARY_FUNCTION(erfc)

TEST_UNARY_FUNCTION(floor)
TEST_UNARY_FUNCTION(round)
TEST_UNARY_FUNCTION(ceil)
TEST_UNARY_FUNCTION(fabs)
TEST_UNARY_FUNCTION(trunc)

TEST_UNARY_FUNCTION(digamma)
TEST(KernelGenerator, log1m_exp_test) {
  MatrixXd m1(3, 3);
  m1 << -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9;

  matrix_cl<double> m1_cl(m1);
  auto tmp = stan::math::log1m_exp(m1_cl);
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = stan::math::log1m_exp(m1);
  EXPECT_NEAR_REL(correct, res);
}
TEST_UNARY_FUNCTION(log1p_exp)
TEST_UNARY_FUNCTION(inv_square)
TEST_UNARY_FUNCTION(inv_logit)
TEST_UNARY_FUNCTION(logit)
TEST_UNARY_FUNCTION(log1m_inv_logit)

#define TEST_CLASSIFICATION_FUNCTION(fun)                                 \
  TEST(KernelGenerator, fun##_test) {                                     \
    using stan::math::fun;                                                \
    MatrixXd m1(3, 3);                                                    \
    m1 << 0.0, 0.2, 0.3, 0.4, 0.5, 0.6, -INFINITY, INFINITY, NAN;         \
                                                                          \
    matrix_cl<double> m1_cl(m1);                                          \
    auto tmp = fun(m1_cl);                                                \
    matrix_cl<bool> res_cl = tmp;                                         \
                                                                          \
    Eigen::Matrix<bool, -1, -1> res = stan::math::from_matrix_cl(res_cl); \
    Eigen::Matrix<bool, -1, -1> correct = fun(m1.array());                \
    EXPECT_TYPED_MATRIX_EQ(correct, res, bool);                           \
  }

TEST_CLASSIFICATION_FUNCTION(isfinite)
TEST_CLASSIFICATION_FUNCTION(isinf)
TEST_CLASSIFICATION_FUNCTION(isnan)

TEST(KernelGenerator, multiple_operations_test) {
  MatrixXd m1(3, 3);
  m1 << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9;

  matrix_cl<double> m1_cl(m1);
  auto tmp = stan::math::exp(stan::math::sin(m1_cl));
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = stan::math::exp(stan::math::sin(m1));
  EXPECT_NEAR_REL(correct, res);
}

TEST(KernelGenerator, multiple_operations_accepts_lvalue_test) {
  MatrixXd m1(3, 3);
  m1 << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9;

  matrix_cl<double> m1_cl(m1);
  auto tmp = stan::math::sin(m1_cl);
  auto tmp2 = stan::math::exp(tmp);
  matrix_cl<double> res_cl = tmp2;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = stan::math::exp(stan::math::sin(m1));
  EXPECT_NEAR_REL(correct, res);
}

TEST(KernelGenerator, multiple_operations_with_includes_test) {
  MatrixXd m1(3, 3);
  m1 << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9;

  matrix_cl<double> m1_cl(m1);
  auto tmp = stan::math::digamma(stan::math::digamma(m1_cl));
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = stan::math::digamma(stan::math::digamma(m1));
  EXPECT_NEAR_REL(correct, res);
}

#endif
