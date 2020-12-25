#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

auto diag_post_multiply_functor
    = [](const auto& a, const auto& b) { return stan::math::diag_post_multiply(a, b); };

TEST(OpenCLMatrixMultiply, prim_rev_values_small) {
  int N = 2;
  int M = 3;

  Eigen::MatrixXd a(N, M);
  a << 1, 2, 3, 4, 5, 6;
  Eigen::VectorXd b(M);
  b << 12, 11, 10;
  stan::math::test::compare_cpu_opencl_prim_rev(diag_post_multiply_functor, a, b);
}

TEST(OpenCLMatrixMultiply, prim_rev_values_N_0) {
  int N = 0;
  int M = 3;

  Eigen::MatrixXd a(N, M);
  Eigen::VectorXd b(M);
  b << 12, 11, 10;
  stan::math::test::compare_cpu_opencl_prim_rev(diag_post_multiply_functor, a, b);
}

TEST(OpenCLMatrixMultiply, prim_rev_values_M_0) {
  int N = 2;
  int M = 0;

  Eigen::MatrixXd a(N, M);
  Eigen::VectorXd b(M);
  stan::math::test::compare_cpu_opencl_prim_rev(diag_post_multiply_functor, a, b);
}

TEST(OpenCLMatrixMultiply, prim_rev_values_large) {
  int N = 71;
  int M = 83;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, M);
  Eigen::VectorXd b = Eigen::VectorXd::Random(M);
  stan::math::test::compare_cpu_opencl_prim_rev(diag_post_multiply_functor, a, b);
}

#endif
