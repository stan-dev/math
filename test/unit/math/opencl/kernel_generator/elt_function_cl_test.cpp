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
TEST_UNARY_FUNCTION(log1m)
TEST_UNARY_FUNCTION(log_inv_logit)

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
TEST_UNARY_FUNCTION(square)
TEST_UNARY_FUNCTION(Phi)
TEST_UNARY_FUNCTION(Phi_approx)
TEST_UNARY_FUNCTION(inv_Phi)

TEST(KernelGenerator, trigamma_test) {
  MatrixXd m1(3, 4);
  m1 << -14.6, -7, -2.2, -1, -0.3, 0, 0.1, 1, 2.2, 7.6, 17.4, 123;

  matrix_cl<double> m1_cl(m1);
  auto tmp = stan::math::trigamma(m1_cl);
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = stan::math::trigamma(m1);
  EXPECT_NEAR_REL(correct, res);
}

#define TEST_CLASSIFICATION_FUNCTION(fun)                                  \
  TEST(KernelGenerator, fun##_test) {                                      \
    using stan::math::fun;                                                 \
    MatrixXd m1(3, 3);                                                     \
    m1 << 0.0, 0.2, 0.3, 0.4, 0.5, 0.6, -INFINITY, INFINITY, NAN;          \
                                                                           \
    matrix_cl<double> m1_cl(m1);                                           \
    auto tmp = fun(m1_cl);                                                 \
    matrix_cl<bool> res_cl = tmp;                                          \
                                                                           \
    Eigen::Matrix<bool, -1, -1> res = stan::math::from_matrix_cl(res_cl);  \
    Eigen::Matrix<bool, -1, -1> correct = fun(m1.array());                 \
    EXPECT_MATRIX_EQ(correct, res);                                        \
                                                                           \
    MatrixXi m1i(3, 3);                                                    \
    m1i << 1, 2, 3, 4, 5, 6, 7, 8, 9;                                      \
    matrix_cl<int> m1i_cl(m1i);                                            \
    auto tmpi = fun(m1_cl);                                                \
    matrix_cl<bool> resi_cl = tmp;                                         \
                                                                           \
    Eigen::Matrix<bool, -1, -1> resi = stan::math::from_matrix_cl(res_cl); \
    Eigen::Matrix<bool, -1, -1> correcti = fun(m1.array());                \
    EXPECT_MATRIX_EQ(correcti, resi);                                      \
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

#define TEST_BINARY_FUNCTION(fun)                                           \
  TEST(KernelGenerator, fun##_test) {                                       \
    MatrixXd m1(3, 3);                                                      \
    m1 << 10.1, 11.2, 12.3, 13.4, 9.5, 24.6, 214.7, 12.8, 3.9;              \
    MatrixXd m2(3, 3);                                                      \
    m2 << 1.1, 2.2, 1.3, 3.4, 4.5, 2.6, 117.2, 11.8, 2.1;                   \
    MatrixXi m3(3, 3);                                                      \
    m3 << 2, 1, 3, 1, 2, 3, 4, 5, 3;                                        \
    MatrixXi m4(3, 3);                                                      \
    m4 << 4, 2, 5, 8, 3, 4, 12, 17, 5;                                      \
    MatrixXi m_size1(3, 2);                                                 \
    m_size1 << 4, 2, 5, 8, 3, 4;                                            \
                                                                            \
    matrix_cl<double> m1_cl(m1);                                            \
    matrix_cl<double> m2_cl(m2);                                            \
    matrix_cl<int> m3_cl(m3);                                               \
    matrix_cl<int> m4_cl(m4);                                               \
    matrix_cl<int> m_size_cl(m_size1);                                      \
                                                                            \
    EXPECT_THROW(stan::math::fun(m1_cl, m_size_cl), std::invalid_argument); \
                                                                            \
    matrix_cl<double> res1_cl = fun(m1_cl, m2_cl);                          \
                                                                            \
    MatrixXd res1 = stan::math::from_matrix_cl(res1_cl);                    \
    MatrixXd correct1 = stan::math::fun(m1, m2);                            \
    EXPECT_NEAR_REL(correct1, res1);                                        \
                                                                            \
    matrix_cl<double> res2_cl = fun(m1_cl, m3_cl);                          \
                                                                            \
    MatrixXd res2 = stan::math::from_matrix_cl(res2_cl);                    \
    MatrixXd correct2 = stan::math::fun(m1, m3);                            \
    EXPECT_NEAR_REL(correct2, res2);                                        \
                                                                            \
    matrix_cl<double> res3_cl = fun(m1_cl, 2.2);                            \
                                                                            \
    MatrixXd res3 = stan::math::from_matrix_cl(res3_cl);                    \
    MatrixXd correct3 = stan::math::fun(m1, 2.2);                           \
    EXPECT_NEAR_REL(correct3, res3);                                        \
                                                                            \
    matrix_cl<double> res4_cl = fun(1000.1, m2_cl);                         \
                                                                            \
    MatrixXd res4 = stan::math::from_matrix_cl(res4_cl);                    \
    MatrixXd correct4 = stan::math::fun(1000.1, m2);                        \
    EXPECT_NEAR_REL(correct4, res4);                                        \
                                                                            \
    matrix_cl<double> res5_cl = fun(m4_cl, m3_cl);                          \
                                                                            \
    MatrixXd res5 = stan::math::from_matrix_cl(res5_cl);                    \
    MatrixXd correct5 = stan::math::fun(m4, m3);                            \
    EXPECT_NEAR_REL(correct5, res5);                                        \
  }

TEST_BINARY_FUNCTION(beta)
TEST_BINARY_FUNCTION(binomial_coefficient_log)
TEST_BINARY_FUNCTION(fdim)
TEST_BINARY_FUNCTION(fmax)
TEST_BINARY_FUNCTION(fmin)
TEST_BINARY_FUNCTION(fmod)
TEST_BINARY_FUNCTION(lbeta)
TEST_BINARY_FUNCTION(lmultiply)
TEST_BINARY_FUNCTION(multiply_log)
TEST_BINARY_FUNCTION(pow)

#endif
