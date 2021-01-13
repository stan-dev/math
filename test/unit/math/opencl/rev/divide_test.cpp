#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

auto divide_functor
    = [](const auto& a, const auto& b) { return stan::math::divide(a, b); };

TEST(OpenCLMatrix_divide, prim_rev_scalar_mat_zero) {
  int M = 0;
  int K = 0;

  double a = 2.0;
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(M, K);
  stan::math::test::compare_cpu_opencl_prim_rev(divide_functor, b, a);
}

TEST(OpenCLMatrix_divide, prim_rev_scalar_mat_values) {
  int N = 71;
  int M = 83;

  double a = 2.0;
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(divide_functor, b, a);
}

#endif
