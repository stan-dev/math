#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <CL/cl2.hpp>
#include <exception>

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using stan::math::matrix_cl;

#define EXPECT_MATRIX_NEAR(A, B, DELTA) \
  for (int i = 0; i < A.size(); i++)    \
    EXPECT_NEAR(A(i), B(i), DELTA);

TEST(MathMatrixCL, operation_cl_errors) {
  EXPECT_THROW(matrix_cl<double> a = stan::math::as_operation_cl(3.5),
               std::domain_error);
  matrix_cl<double> b(3, 3);
  matrix_cl<double> c(4, 3);
  EXPECT_THROW((b + 3).evaluate_into(c), std::invalid_argument);
}

TEST(MathMatrixCL, kernel_caching) {
  MatrixXd m1(3, 3);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  MatrixXd m2(3, 3);
  m2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);
  auto tmp = m1_cl + 0.1234 * m2_cl;
  EXPECT_EQ(decltype(tmp)::cache<matrix_cl<double>>::kernel(), nullptr);
  matrix_cl<double> res_cl = tmp;
  cl_kernel cached_kernel = decltype(tmp)::cache<matrix_cl<double>>::kernel();
  EXPECT_NE(cached_kernel, nullptr);

  auto tmp2 = m1_cl + 0.1234 * m2_cl;
  matrix_cl<double> res2_cl = tmp2;
  EXPECT_EQ(decltype(tmp)::cache<matrix_cl<double>>::kernel(), cached_kernel);

  matrix_cl<double> res3_cl = res_cl + 0.1234 * res2_cl;
  EXPECT_EQ(decltype(tmp)::cache<matrix_cl<double>>::kernel(), cached_kernel);

  EXPECT_EQ(decltype(tmp)::cache<matrix_cl<int>>::kernel(), nullptr);
}

#endif
