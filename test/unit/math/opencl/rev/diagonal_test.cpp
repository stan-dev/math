#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto diagonal_functor
    = [](const auto& a) { return stan::math::diagonal(a).eval(); };

TEST(OpenCLMatrixMultiply, prim_rev_values_small) {
  int N = 2;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(diagonal_functor, a);
}

TEST(OpenCLMatrixMultiply, prim_rev_values_N_0) {
  int N = 0;

  Eigen::MatrixXd a(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(diagonal_functor, a);
}

TEST(OpenCLMatrixMultiply, prim_rev_values_large) {
  int N = 71;
  
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(diagonal_functor, a);
}

#endif
