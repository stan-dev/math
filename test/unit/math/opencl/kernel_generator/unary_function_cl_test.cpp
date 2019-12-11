#ifdef STAN_OPENCL

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/opencl/kernel_generator/unary_function_cl.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <test/unit/math/opencl/kernel_generator/reference_kernel.hpp>
#include <stan/math.hpp>
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

#define EXPECT_MATRIX_NEAR(A, B, DELTA) \
  for (int i = 0; i < A.size(); i++)    \
    EXPECT_NEAR(A(i), B(i), DELTA);

#define TEST_FUNCTION(fun)                             \
  TEST(MathMatrixCL, fun##_test) {                     \
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
    EXPECT_MATRIX_NEAR(correct, res, 1e-9);            \
  }

TEST_FUNCTION(rsqrt)
TEST_FUNCTION(sqrt)
TEST_FUNCTION(cbrt)

TEST_FUNCTION(exp)
TEST_FUNCTION(exp2)
TEST_FUNCTION(expm1)
TEST_FUNCTION(log)
TEST_FUNCTION(log2)
TEST_FUNCTION(log10)
TEST_FUNCTION(log1p)

TEST_FUNCTION(sin)
TEST_FUNCTION(sinh)
TEST_FUNCTION(cos)
TEST_FUNCTION(cosh)
TEST_FUNCTION(tan)
TEST_FUNCTION(tanh)

TEST_FUNCTION(asin)
TEST_FUNCTION(asinh)
TEST_FUNCTION(acos)

TEST(MathMatrixCL, acosh_test) {
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
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

TEST_FUNCTION(atan)
TEST_FUNCTION(atanh)

TEST_FUNCTION(tgamma)
TEST_FUNCTION(lgamma)
TEST_FUNCTION(erf)
TEST_FUNCTION(erfc)

TEST_FUNCTION(floor)
TEST_FUNCTION(round)
TEST_FUNCTION(ceil)

TEST(MathMatrixCL, multiple_operations_test) {
  MatrixXd m1(3, 3);
  m1 << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9;

  matrix_cl<double> m1_cl(m1);
  auto tmp = stan::math::exp(stan::math::sin(m1_cl));
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = stan::math::exp(stan::math::sin(m1));
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

TEST(MathMatrixCL, multiple_operations_accepts_lvalue_test) {
  MatrixXd m1(3, 3);
  m1 << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9;

  matrix_cl<double> m1_cl(m1);
  auto tmp = stan::math::sin(m1_cl);
  auto tmp2 = stan::math::exp(tmp);
  matrix_cl<double> res_cl = tmp2;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = stan::math::exp(stan::math::sin(m1));
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}
#endif
