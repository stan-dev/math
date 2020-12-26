#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

auto diag_pre_multiply_functor = [](const auto& a, const auto& b) {
  return stan::math::diag_pre_multiply(a, b);
};

TEST(OpenCL_diag_pre_multiply, prim_rev_values_small) {
  int N = 2;
  int M = 3;

  Eigen::VectorXd a(N);
  a << 12, 11;
  Eigen::MatrixXd b(N, M);
  b << 1, 2, 3, 4, 5, 6;
  stan::math::test::compare_cpu_opencl_prim_rev(diag_pre_multiply_functor, a,
                                                b);
}

TEST(OpenCL_diag_pre_multiply, prim_rev_values_N_0) {
  int N = 0;
  int M = 3;

  Eigen::VectorXd a(N);
  Eigen::MatrixXd b(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(diag_pre_multiply_functor, a,
                                                b);
}

TEST(OpenCL_diag_pre_multiply, prim_rev_values_M_0) {
  int N = 2;
  int M = 0;

  Eigen::VectorXd a(N);
  a << 12, 11;
  Eigen::MatrixXd b(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(diag_pre_multiply_functor, a,
                                                b);
}

TEST(OpenCL_diag_pre_multiply, prim_rev_values_large) {
  int N = 71;
  int M = 83;

  Eigen::VectorXd a = Eigen::VectorXd::Random(N);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(diag_pre_multiply_functor, a,
                                                b);
}

#endif
