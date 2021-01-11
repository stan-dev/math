#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto Phi_approx_functor
    = [](const auto& a) { return stan::math::Phi_approx(a); };

TEST(OpenCL_Phi_approx, prim_rev_values_small) {
  Eigen::VectorXd a(8);
  a << -0.22, -0.8, 0.5, 1, 0.15, 0.3, 0.34, -0.04;
  stan::math::test::compare_cpu_opencl_prim_rev(Phi_approx_functor, a);
}

TEST(OpenCL_Phi_approx, prim_rev_size_0) {
  int N = 0;

  Eigen::MatrixXd a(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(Phi_approx_functor, a);
}

TEST(OpenCL_Phi_approx, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(Phi_approx_functor, a);
}

#endif
