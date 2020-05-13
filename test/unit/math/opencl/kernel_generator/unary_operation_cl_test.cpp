#ifdef STAN_OPENCL

#include <stan/math/prim.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <test/unit/math/opencl/kernel_generator/reference_kernel.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <string>

#define EXPECT_MATRIX_NEAR(A, B, DELTA) \
  for (int i = 0; i < A.size(); i++)    \
    EXPECT_NEAR(A(i), B(i), DELTA);

TEST(KernelGenerator, logical_negation_test) {
  using stan::math::matrix_cl;
  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> m1(3, 3);
  m1 << true, false, true, true, false, false, true, false, false;

  matrix_cl<bool> m1_cl(m1);
  matrix_cl<bool> res_cl = !m1_cl;

  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> res
      = stan::math::from_matrix_cl(res_cl);
  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> correct = !m1.array();
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, unary_minus_test) {
  using stan::math::matrix_cl;
  Eigen::MatrixXd m1(3, 3);
  m1 << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> res_cl = -m1_cl;

  Eigen::MatrixXd res = stan::math::from_matrix_cl(res_cl);
  Eigen::MatrixXd correct = -m1;
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

#endif
