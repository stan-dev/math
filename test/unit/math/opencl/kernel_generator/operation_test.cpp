#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/zeros.hpp>
#include <stan/math/opencl/copy.hpp>
#include <test/unit/util.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <CL/cl2.hpp>
#include <exception>

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using stan::math::matrix_cl;

TEST(KernelGenerator, operation_cl_errors) {
  EXPECT_THROW(matrix_cl<double> a = stan::math::as_operation_cl(3.5),
               std::domain_error);
  EXPECT_THROW(stan::math::as_operation_cl(3.5).eval(), std::domain_error);
}

TEST(KernelGenerator, kernel_caching) {
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

TEST(MathMatrixCL, events_write_after_write) {
  using stan::math::matrix_cl;
  matrix_cl<double> zero_cl(3, 3);
  zero_cl.zeros();
  zero_cl.wait_for_read_write_events();

  for (int j = 0; j < 3000; j++) {
    matrix_cl<double> m_cl(3, 3);

    for (int i = 0; i < 4; i++) {
      m_cl = zero_cl + i;
    }

    Eigen::MatrixXd res = stan::math::from_matrix_cl(m_cl);
    Eigen::MatrixXd correct = Eigen::MatrixXd::Constant(3, 3, 3);

    EXPECT_MATRIX_NEAR(res, correct, 1e-13);
  }
}

TEST(MathMatrixCL, events_read_after_write_and_write_after_read) {
  using stan::math::matrix_cl;
  int iters = 3000;

  matrix_cl<double> m1_cl(3, 3);
  matrix_cl<double> m2_cl(3, 3);
  m1_cl.zeros();

  for (int j = 0; j < iters; j++) {
    m2_cl = m1_cl + 1;
    m1_cl = m2_cl + 1;
  }
  Eigen::MatrixXd res = stan::math::from_matrix_cl(m1_cl);
  Eigen::MatrixXd correct = Eigen::MatrixXd::Constant(3, 3, 2 * iters);

  EXPECT_MATRIX_NEAR(res, correct, 1e-13);
}

#endif
