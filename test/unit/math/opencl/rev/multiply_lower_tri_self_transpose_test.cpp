#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto multiply_lower_tri_self_transpose_functor = [](const auto& a) {
  return stan::math::multiply_lower_tri_self_transpose(a);
};

TEST(OpenCLMultiplyLowerTriSelfTranspose, prim_rev_values_small) {
  Eigen::MatrixXd a(4, 2);
  a << -2.2, -0.8, 0.5, 1, 1.5, 3, 3.4, 4;
  stan::math::test::compare_cpu_opencl_prim_rev(
      multiply_lower_tri_self_transpose_functor, a);
}

TEST(OpenCLMultiplyLowerTriSelfTranspose, prim_rev_size_0) {
  int N = 0;

  Eigen::MatrixXd a(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multiply_lower_tri_self_transpose_functor, a);
}

TEST(OpenCLMultiplyLowerTriSelfTranspose, prim_rev_values_large) {
  int N = 71;
  int M = 87;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(
      multiply_lower_tri_self_transpose_functor, a);
}

#endif
