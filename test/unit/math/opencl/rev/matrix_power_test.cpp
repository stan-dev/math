#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto matrix_power_functor
    = [](const auto& a, int b) { return stan::math::matrix_power(a, b); };

TEST(OpenCLMatrixPower, prim_rev_values_small) {
  Eigen::MatrixXd a(3, 3);
  a << -2.2, -0.8, 0.5, 1, 1.5, 3, 3.4, 4, 0.4;
  stan::math::test::compare_cpu_opencl_prim_rev(matrix_power_functor, a, 7);
}

TEST(OpenCLMatrixPower, prim_rev_size_0) {
  int N = 0;

  Eigen::MatrixXd a(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(matrix_power_functor, a, 3);
}

TEST(OpenCLMatrixPower, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(matrix_power_functor, a, 17);
}

#endif
