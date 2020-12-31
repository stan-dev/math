#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto erfc_functor = [](const auto& a) { return stan::math::erfc(a); };

TEST(OpenCLerfc, prim_rev_values_small) {
  Eigen::VectorXd a(7);
  a << -2.6, -2, -1, -0.2, 1, 1.3, 2.6;
  stan::math::test::compare_cpu_opencl_prim_rev(erfc_functor, a);
}

TEST(OpenCLerfc, prim_rev_size_0) {
  int N = 0;

  Eigen::MatrixXd a(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(erfc_functor, a);
}

TEST(OpenCLerfc, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(erfc_functor, a);
}

#endif
