#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <gtest/gtest.h>

auto transpose_functor
    = [](const auto& a) { return stan::math::transpose(a).eval(); };

TEST(OpenCLMatrixTranspose, prim_rev_values_small) {
  int N = 2;
  int M = 3;

  Eigen::MatrixXd a(N, M);
  a << 1, 2, 3, 4, 5, 6;
  stan::math::test::compare_cpu_opencl_prim_rev(transpose_functor, a);
}

TEST(OpenCLMatrixTranspose, prim_rev_values_N_0) {
  int N = 0;
  int M = 2;

  Eigen::MatrixXd a(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(transpose_functor, a);
}

TEST(OpenCLMatrixTranspose, prim_rev_values_M_0) {
  int N = 2;
  int M = 0;

  Eigen::MatrixXd a(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(transpose_functor, a);
}

TEST(OpenCLMatrixTranspose, prim_rev_values_large) {
  int N = 71;
  int M = 83;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(transpose_functor, a);
}

#endif
