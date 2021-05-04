#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

auto log_diff_exp_functor = [](const auto& a, const auto& b) {
  return stan::math::log_diff_exp(a, b);
};

TEST(OpenCLMatrix_log_diff_exp, prim_rev_values_small) {
  int N = 2;
  int M = 3;

  Eigen::MatrixXd a(N, M);
  a << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd b(N, M);
  b << 12, 11, 10, 9, 8, 7;
  stan::math::test::compare_cpu_opencl_prim_rev(log_diff_exp_functor, a, b);
}

TEST(OpenCLMatrix_log_diff_exp, prim_rev_values_M_0) {
  int N = 2;
  int M = 0;

  Eigen::MatrixXd a(N, M);
  Eigen::MatrixXd b(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(log_diff_exp_functor, a, b);

  Eigen::MatrixXd c(M, N);
  Eigen::MatrixXd d(M, N);
  stan::math::test::compare_cpu_opencl_prim_rev(log_diff_exp_functor, c, d);
}

TEST(OpenCLMatrix_log_diff_exp, prim_rev_values_large) {
  int N = 71;
  int M = 83;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, M);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(log_diff_exp_functor, a, b);
}

TEST(OpenCLMatrix_log_diff_exp, prim_rev_scalar_values_large) {
  int N = 71;
  int M = 83;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, M);
  double b = 0.3;
  stan::math::test::compare_cpu_opencl_prim_rev(log_diff_exp_functor, a, b);
  stan::math::test::compare_cpu_opencl_prim_rev(log_diff_exp_functor, b, a);
}

#endif
