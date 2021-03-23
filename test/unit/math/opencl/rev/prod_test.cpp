#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto prod_functor = [](const auto& a) { return stan::math::prod(a); };

TEST(OpenCLProd, prim_rev_values_small) {
  int N = 2;
  int M = 3;

  Eigen::MatrixXd a(N, M);
  a << 1, 2, 3, 4, 5, 6;
  stan::math::test::compare_cpu_opencl_prim_rev(prod_functor, a);
}

TEST(OpenCLProd, prim_rev_values_zero_rows) {
  int N = 0;
  int M = 3;

  Eigen::MatrixXd a(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(prod_functor, a);
}

TEST(OpenCLProd, prim_rev_values_zero_cols) {
  int N = 2;
  int M = 0;

  Eigen::MatrixXd a(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(prod_functor, a);
}

TEST(OpenCLProd, prim_rev_values_large) {
  int N = 71;
  int M = 83;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(prod_functor, a);
}

#endif
