#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto inv_Phi_functor = [](const auto& a) { return stan::math::inv_Phi(a); };

TEST(OpenCL_inv_Phi, prim_rev_values_small) {
  Eigen::VectorXd a(7);
  a << 0.02425, 0.97575, 0.01, 0.1, 0.98, 0.5, 1.0;
  stan::math::test::compare_cpu_opencl_prim_rev(inv_Phi_functor, a);
}

TEST(OpenCL_inv_Phi, prim_rev_size_0) {
  int N = 0;

  Eigen::MatrixXd a(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(inv_Phi_functor, a);
}

TEST(OpenCL_inv_Phi, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N).array().abs();
  stan::math::test::compare_cpu_opencl_prim_rev(inv_Phi_functor, a);
}

#endif
