#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto tgamma_functor = [](const auto& a) { return stan::math::tgamma(a); };

TEST_F(OpenCLRevTests, _tgamma_prim_rev_values_small) {
  Eigen::VectorXd a(14);
  a << -15.2, -10.1, -0.5, 0.5, 1, 1.0, 1.3, 5, 10, -2.6, -2.9, -0.2, 1.3, 3;
  stan::math::test::compare_cpu_opencl_prim_rev(tgamma_functor, a);
}

TEST_F(OpenCLRevTests, _tgamma_prim_rev_size_0) {
  int N = 0;

  Eigen::MatrixXd a(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(tgamma_functor, a);
}

TEST_F(OpenCLRevTests, _tgamma_prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(tgamma_functor, a);
}

#endif
