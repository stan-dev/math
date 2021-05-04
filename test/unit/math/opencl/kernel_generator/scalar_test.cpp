#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using stan::math::matrix_cl;

TEST(KernelGenerator, scalar_test) {
  MatrixXd m1(3, 3);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;

  matrix_cl<double> m1_cl(m1);
  auto tmp = m1_cl + 0.1234;
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  EXPECT_EQ(1.1234, res(0, 0));
  EXPECT_EQ(2.1234, res(0, 1));
  EXPECT_EQ(3.1234, res(0, 2));
  EXPECT_EQ(4.1234, res(1, 0));
  EXPECT_EQ(5.1234, res(1, 1));
  EXPECT_EQ(6.1234, res(1, 2));
  EXPECT_EQ(7.1234, res(2, 0));
  EXPECT_EQ(8.1234, res(2, 1));
  EXPECT_EQ(9.1234, res(2, 2));
}

TEST(KernelGenerator, scalar_test_accept_lvalue) {
  MatrixXd m1(3, 3);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;

  matrix_cl<double> m1_cl(m1);
  double scal = 0.1234;
  auto tmp = m1_cl + scal;
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  EXPECT_EQ(1.1234, res(0, 0));
  EXPECT_EQ(2.1234, res(0, 1));
  EXPECT_EQ(3.1234, res(0, 2));
  EXPECT_EQ(4.1234, res(1, 0));
  EXPECT_EQ(5.1234, res(1, 1));
  EXPECT_EQ(6.1234, res(1, 2));
  EXPECT_EQ(7.1234, res(2, 0));
  EXPECT_EQ(8.1234, res(2, 1));
  EXPECT_EQ(9.1234, res(2, 2));
}

#endif
