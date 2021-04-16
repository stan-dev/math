#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <test/unit/math/opencl/kernel_generator/reference_kernel.hpp>
#include <test/unit/util.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <string>

using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using stan::math::matrix_cl;

TEST(KernelGenerator, cast_zero_size) {
  matrix_cl<double> m_zero(0, 0);

  EXPECT_NO_THROW(stan::math::cast<int>(m_zero));
}

TEST(KernelGenerator, cast_test) {
  using stan::math::cast;
  MatrixXd m1(2, 3);
  m1 << 1, 2.5, 3.4, 4.7, 5.9, 6.3;

  matrix_cl<double> m1_cl(m1);

  matrix_cl<int> res_cl = cast<int>(m1_cl);

  MatrixXi res = stan::math::from_matrix_cl(res_cl);

  MatrixXi correct = m1.cast<int>();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, cast_multiple_operations) {
  using stan::math::cast;
  MatrixXd m1(2, 3);
  m1 << 1, 2.5, 3.4, 4.7, 5.9, 6.3;

  matrix_cl<double> m1_cl(m1);
  auto tmp = cast<double>(cast<int>(m1_cl));
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m1.cast<int>().template cast<double>();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, cast_multiple_operations_lvalue) {
  using stan::math::cast;
  MatrixXd m1(2, 3);
  m1 << 1, 2.5, 3.4, 4.7, 5.9, 6.3;

  matrix_cl<double> m1_cl(m1);
  auto tmp = cast<int>(m1_cl);
  matrix_cl<double> res_cl = cast<double>(tmp);

  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m1.cast<int>().template cast<double>();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

#endif
