#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto lgamma_functor = [](const auto& a) { return stan::math::lgamma(a); };

TEST(OpenCL_lgamma, prim_rev_values_small) {
  Eigen::VectorXd a(9);
  a << -15.2, -10.2, -0.5, 0.5, 1, 1.0, 1.3, 5, 10;
  stan::math::test::compare_cpu_opencl_prim_rev(lgamma_functor, a);
}

TEST(OpenCL_lgamma, prim_rev_size_0) {
  int N = 0;

  Eigen::MatrixXd a(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(lgamma_functor, a);
}

TEST(OpenCL_lgamma, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a
      = Eigen::MatrixXd::Constant(N, N, 1.1) + Eigen::MatrixXd::Random(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(lgamma_functor, a);
}

#endif
