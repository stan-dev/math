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
}

TEST(MathMatrixCL, kernel_caching) {
  MatrixXd m1(3, 3);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  MatrixXd m2(3, 3);
  m2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);
  auto tmp = m1_cl + 0.1234 * m2_cl;
  using cache = stan::math::internal::multi_result_kernel_internal<
      0, stan::math::load_<matrix_cl<double>&>>::inner<const decltype(tmp)&>;
  using unused_cache = stan::math::internal::multi_result_kernel_internal<
      0, stan::math::load_<matrix_cl<int>&>>::inner<const decltype(tmp)&>;
  EXPECT_EQ(cache::kernel_(), nullptr);
  matrix_cl<double> res_cl = tmp;
  cl_kernel cached_kernel = cache::kernel_();
  EXPECT_NE(cached_kernel, nullptr);

  auto tmp2 = m1_cl + 0.1234 * m2_cl;
  matrix_cl<double> res2_cl = tmp2;
  EXPECT_EQ(cache::kernel_(), cached_kernel);

  matrix_cl<double> res3_cl = res_cl + 0.1234 * res2_cl;
  EXPECT_EQ(cache::kernel_(), cached_kernel);

  EXPECT_EQ(unused_cache::kernel_(), nullptr);
}

#endif
