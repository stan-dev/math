#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto sd_functor = [](const auto& a) { return stan::math::sd(a); };

TEST(OpenCLSd, prim_rev_values_small) {
  Eigen::VectorXd a(8);
  a << -2.2, -0.8, 0.5, 1, 1.5, 3, 3.4, 4;
  stan::math::test::compare_cpu_opencl_prim_rev(sd_functor, a);
}

TEST(OpenCLSd, prim_rev_size_1) {
  int N = 1;

  Eigen::MatrixXd a(N, N);
  a << 123;
  stan::math::test::compare_cpu_opencl_prim_rev(sd_functor, a);
}

TEST(OpenCLSd, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(sd_functor, a);
}

#endif
