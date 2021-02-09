#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>

auto beta_functor
    = [](const auto& a, const auto& b) { return stan::math::beta(a, b); };

TEST(OpenCL_beta, beta_small) {
  Eigen::VectorXd in1(4);
  in1 << 0.5, 3.4, 5.2, 7.5;
  Eigen::VectorXd in2(4);
  in2 << 3.3, 0.9, 6.7, 1.8;
  stan::math::test::compare_cpu_opencl_prim_rev(beta_functor, in1, in2);
}

TEST(OpenCL_beta, zero) {
  Eigen::VectorXd in1;
  Eigen::VectorXd in2;
  stan::math::test::compare_cpu_opencl_prim_rev(beta_functor, in1, in2);
}

TEST(OpenCL_beta, prim_rev_values_large) {
  int N = 500;

  Eigen::VectorXd a = Eigen::VectorXd::Random(N);
  Eigen::VectorXd b = Eigen::VectorXd::Random(N);
  stan::math::test::compare_cpu_opencl_prim_rev(beta_functor, a, b);
}

TEST(OpenCL_beta, prim_rev_scalar_values_large) {
  int N = 71;
  int M = 83;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, M);
  double b = 0.3;
  stan::math::test::compare_cpu_opencl_prim_rev(beta_functor, a, b);
  stan::math::test::compare_cpu_opencl_prim_rev(beta_functor, b, a);
}
#endif
