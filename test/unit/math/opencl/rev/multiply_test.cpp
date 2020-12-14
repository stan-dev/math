#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

auto matrix_multiply_functor
    = [](const auto& a, const auto& b) { return stan::math::multiply(a, b); };

TEST(OpenCLMatrixMultiply, prim_rev_values_small) {
  int N = 2;
  int M = 3;
  int K = 4;

  Eigen::MatrixXd a(N, M);
  a << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd b(M, K);
  b << 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1;
  stan::math::test::compare_cpu_opencl_prim_rev(matrix_multiply_functor, a, b);
}

TEST(OpenCLMatrixMultiply, prim_rev_values_N_0) {
  int N = 0;
  int M = 3;
  int K = 4;

  Eigen::MatrixXd a(N, M);
  Eigen::MatrixXd b(M, K);
  b << 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1;
  stan::math::test::compare_cpu_opencl_prim_rev(matrix_multiply_functor, a, b);
}

TEST(OpenCLMatrixMultiply, prim_rev_values_M_0) {
  int N = 2;
  int M = 0;
  int K = 4;

  Eigen::MatrixXd a(N, M);
  Eigen::MatrixXd b(M, K);
  stan::math::test::compare_cpu_opencl_prim_rev(matrix_multiply_functor, a, b);
}

TEST(OpenCLMatrixMultiply, prim_rev_values_K_0) {
  int N = 2;
  int M = 3;
  int K = 0;

  Eigen::MatrixXd a(N, M);
  a << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd b(M, K);
  stan::math::test::compare_cpu_opencl_prim_rev(matrix_multiply_functor, a, b);
}

TEST(OpenCLMatrixMultiply, prim_rev_values_large) {
  int N = 71;
  int M = 83;
  int K = 97;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, M);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(M, K);
  stan::math::test::compare_cpu_opencl_prim_rev(matrix_multiply_functor, a, b);
}

TEST(OpenCLMatrixMultiply, prim_rev_scalar_mat_zero) {
  int M = 0;
  int K = 0;

  double a = 2.0;
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(M, K);
  stan::math::test::compare_cpu_opencl_prim_rev(matrix_multiply_functor, a, b);
  stan::math::test::compare_cpu_opencl_prim_rev(matrix_multiply_functor, b, a);
}

TEST(OpenCLMatrixMultiply, prim_rev_scalar_mat_values) {
  int N = 71;
  int M = 83;

  double a = 2.0;
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(matrix_multiply_functor, a, b);
  stan::math::test::compare_cpu_opencl_prim_rev(matrix_multiply_functor, b, a);
}

#endif
