#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto cholesky_decompose_functor
    = [](const auto& a) { return stan::math::cholesky_decompose(a); };

TEST(OpenCLCholeskyDecompose, prim_rev_values_small) {
  int N = 3;

  Eigen::MatrixXd a(N, N);
  a << 9, 2, 3, 2, 7, 5, 3, 5, 10;
  stan::math::test::compare_cpu_opencl_prim_rev(cholesky_decompose_functor, a);
}

TEST(OpenCLCholeskyDecompose, prim_rev_size_1) {
  int N = 1;

  Eigen::MatrixXd a(N, N);
  a << 9;
  stan::math::test::compare_cpu_opencl_prim_rev(cholesky_decompose_functor, a);
}

TEST(OpenCLCholeskyDecompose, prim_rev_size_0) {
  int N = 0;

  Eigen::MatrixXd a(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(cholesky_decompose_functor, a);
}

TEST(OpenCLCholeskyDecompose, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N);
  a = a + a.transpose().eval();
  a.diagonal().array() += 2 * N;
  stan::math::test::compare_cpu_opencl_prim_rev(cholesky_decompose_functor, a);

  N = 1251;

  Eigen::MatrixXd b = Eigen::MatrixXd::Random(N, N);
  b = b + b.transpose().eval();
  b.diagonal().array() += 2 * N;
  stan::math::test::compare_cpu_opencl_prim_rev(cholesky_decompose_functor, b);

  N = 1704;

  Eigen::MatrixXd c = Eigen::MatrixXd::Random(N, N);
  c = c + c.transpose().eval();
  c.diagonal().array() += 2 * N;
  stan::math::test::compare_cpu_opencl_prim_rev(cholesky_decompose_functor, c);
}

#endif
