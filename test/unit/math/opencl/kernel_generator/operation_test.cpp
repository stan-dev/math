#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator/binary_operation.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using stan::math::matrix_cl;


#define EXPECT_MATRIX_NEAR(A, B, DELTA) \
  for (int i = 0; i < A.size(); i++)    \
    EXPECT_NEAR(A(i), B(i), DELTA);

TEST(MathMatrixCL, reuse_expression){
  MatrixXd m1(3, 3);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  MatrixXd m2(3, 3);
  m2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);
  auto tmp = m1_cl + 0.1234 + m2_cl;
  matrix_cl<double> res_cl = stan::math::elewise_multiplication(tmp - 6, tmp);

  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  auto tmp_eig = m1.array() + 0.1234 + m2.array();
  MatrixXd res_eig = (tmp_eig - 6)*tmp_eig;

  EXPECT_MATRIX_NEAR(res_eig,res,1e-9);
}

TEST(MathMatrixCL, reuse_kernel){
  MatrixXd m1(3, 3);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  MatrixXd m2(3, 3);
  m2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);
  auto tmp = m1_cl + 0.1234 * m2_cl;
  matrix_cl<double> res_cl = tmp;

  auto tmp2 = m1_cl + 0.1234 * m2_cl;
  matrix_cl<double> res2_cl = tmp2;

  matrix_cl<double> res3_cl = res_cl + 0.1234 * res2_cl;
}

#endif
