#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto to_array_2d_functor = [](const auto& a) {
  return stan::math::to_matrix(stan::math::to_array_2d(a));
};

TEST(OpenCLToArray2D, prim_rev_values_small) {
  Eigen::MatrixXd a(3, 2);
  stan::math::test::compare_cpu_opencl_prim_rev(to_array_2d_functor, a);
}

TEST(OpenCLToArray2D, prim_rev_size_0) {
  Eigen::MatrixXd a(0, 0);
  stan::math::test::compare_cpu_opencl_prim_rev(to_array_2d_functor, a);
}

TEST(OpenCLToArray2D, prim_rev_values_large) {
  int N = 71;
  int M = 87;
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(to_array_2d_functor, a);
}

#endif
