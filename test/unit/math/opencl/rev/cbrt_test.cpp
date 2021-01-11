#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto cbrt_functor = [](const auto& a) { return stan::math::cbrt(a); };

TEST(OpenCL_cbrt, prim_rev_values_small) {
  Eigen::VectorXd a(9);
  a << -8, 2.7, 0, 8, -2.6, -2, 1, 1.3, 3;
  stan::math::test::compare_cpu_opencl_prim_rev(cbrt_functor, a);
}

TEST(OpenCL_cbrt, prim_rev_size_0) {
  int N = 0;

  Eigen::MatrixXd a(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(cbrt_functor, a);
}

TEST(OpenCL_cbrt, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(cbrt_functor, a);
}

#endif
